version 1.0

# Made by Fabio Cunial. All references to "hgsvc2" are purely for the purpose of maximizing confusion.


# Performs a trivial bcftools merge of multiple samples with multiple callers
# per sample. The output VCF contains no multiallelic records and no exact
# duplicates.
#
# Remark: task $InterSampleMerge$ should be made parallel by chromosome.
#
workflow BcftoolsMergeNoDuplicates {
    input {
        Array[File] vcf_gz
        Array[File] tbi
        File reference_fa
        File reference_fai
        Int ram_gb_intersample = 200
        Int compression_level = 1
    }
    parameter_meta {
    }

    scatter(i in range(length(vcf_gz))) {
        call IntraSampleMerge {
            input:
                sample_vcf_gz = [vcf_gz[i]],
                sample_tbi = [tbi[i]],
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                compression_level = compression_level
        }
    }
    call InterSampleMerge {
        input:
            input_vcf_gz = IntraSampleMerge.output_vcf_gz,
            input_tbi = IntraSampleMerge.output_tbi,
            ram_gb = ram_gb_intersample,
            compression_level = compression_level
    }

    output {
        File output_vcf_gz = InterSampleMerge.output_vcf_gz
        File output_tbi = InterSampleMerge.output_tbi
    }
}


# Remark: we use $bcftools merge$ instead of $bcftools concat$, since we must
# collapse identical calls made by different callers (otherwise they would
# remain in the inter-sample VCF, since the inter-sample $bcftools merge$ does
# not collapse records from the same sample).
#
# Remark: $bcftools merge$ collapses into the same record every record with the
# same $CHR,POS,REF,ALT$, disregarding the INFO field and in particular 
# differences in SVLEN and END. This may delete information for symbolic
# ALTs. Our script makes sure that only symbolic records with the same SVLEN and
# END are collapsed into the same record.
#
# Remark: we do not consider STRAND in the above.
#
# Remark: symbolic ALTs in the output VCF are not necessarily identical to the
# symbolic ALTs in the input, and the input genotypes are discarded (i.e. the
# output contains artificial FORMAT and SAMPLE columns where all calls are 0/1).
#
# Performance on each AoU 8x sample:
# COMMAND           RUNTIME     N_CPUS      MAX_RSS
# bcftools sort     4m          1           230M
# bcftools norm     60s         2           230M
# awk               30s         1           156M
# rm-dup            2s          1           16M
# bcftools merge    3m          2           250M
# bcftools norm     3m          2           230M
# bgzip             50s         2           9M
#
task IntraSampleMerge {
    input {
        Array[File] sample_vcf_gz
        Array[File] sample_tbi
        File reference_fa
        File reference_fai
        Int compression_level = 1
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(sample_vcf_gz, "GB")) + 10
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    Int ram_gb = 4
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_gb} - 2 ))
        
        function cleanVCF() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_VCF=$2
            
            # - Ensuring that the input file is sorted
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${INPUT_VCF_GZ} > tmp0.vcf.gz
            tabix -f tmp0.vcf.gz
            
            # - Removing multiallelic records.
            # - Fixing wrong REF values (which may occur e.g. in sniffles).
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z tmp0.vcf.gz > tmp1.vcf.gz
            tabix -f tmp1.vcf.gz
            rm -f tmp0.vcf.gz*
            
            # - Storing SVLEN and END in symbolic ALTs
            bcftools view --header-only tmp1.vcf.gz > tmp2.vcf
            ${TIME_COMMAND} bcftools view --no-header tmp1.vcf.gz | awk '{ \
                tag="artificial"; \
                if ($5=="<DEL>" || $5=="<INS>" || $5=="<INV>" || $5=="<DUP>" || $5=="<CNV>") { \
                    svtype=substr($5,2,3); \
                    end=""; \
                    svlen=""; \
                    n=split($8,A,";"); \
                    for (i=1; i<=n; i++) { \
                        if (substr(A[i],1,4)=="END=") end=substr(A[i],5); \
                        else if (substr(A[i],1,6)=="SVLEN=") { \
                            if (substr(A[i],7,1)=="-") svlen=substr(A[i],8); \
                            else svlen=substr(A[i],7); \
                        } \
                    } \
                    $5="<" svtype ":" tag ":" (length(end)==0?"?":end) ":" (length(svlen)==0?"?":svlen) ">" \
                }; \
                printf("%s",$1); \
                for (i=2; i<=NF; i++) printf("\t%s",$i); \
                printf("\n"); \
            }' >> tmp2.vcf
            rm -f tmp1.vcf.gz*
            ${TIME_COMMAND} bgzip --threads ${N_THREADS} --compress-level ~{compression_level} tmp2.vcf
            tabix -f tmp2.vcf.gz
            
            # - Removing identical records
            # See <https://github.com/samtools/bcftools/issues/1089>.
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --rm-dup exact --output-type v tmp2.vcf.gz > ${OUTPUT_VCF}
            rm -f tmp2.vcf.gz*
            ${TIME_COMMAND} bgzip --threads ${N_THREADS} --compress-level ~{compression_level} ${OUTPUT_VCF}
            tabix -f ${OUTPUT_VCF}.gz
        }
        
        # Merging all single-caller VCFs
        INPUT_FILES=~{sep=',' sample_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        i="0"; SAMPLE_ID="";
        for INPUT_FILE in ${INPUT_FILES}; do
            SAMPLE_ID=$(bcftools view --header-only ${INPUT_FILE} | tail -n 1 | cut -f 10)
            cleanVCF ${INPUT_FILE} ${i}.vcf
            echo ${i}.vcf.gz >> list.txt
            rm -f ${INPUT_FILE}
            i=$(( ${i} + 1 ))
        done
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --force_single --file-list list.txt --output-type z > tmp2.vcf.gz
        tabix -f tmp2.vcf.gz
        
        # Removing multiallelic records, if any are generated during the merge.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z tmp2.vcf.gz > tmp3.vcf.gz
        tabix -f tmp3.vcf.gz
        rm -f tmp2.vcf.gz*
        
        # Removing GTs.
        bcftools view --header-only tmp3.vcf.gz > header.txt
        N_ROWS=$(wc -l < header.txt)
        head -n $(( ${N_ROWS} - 1 )) header.txt > out.vcf
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${SAMPLE_ID}" >> out.vcf
        ${TIME_COMMAND} bcftools view --no-header tmp3.vcf.gz | awk '{ \
            printf("%s",$1); \
            for (i=2; i<=8; i++) printf("\t%s",$i); \
            printf("\tGT\t0/1\n"); \
        }' >> out.vcf
        ${TIME_COMMAND} bgzip --threads ${N_THREADS} --compress-level ~{compression_level} out.vcf
        tabix -f out.vcf.gz
        ls -laht; tree
    >>>

    output {
        File output_vcf_gz = work_dir + "/out.vcf.gz"
        File output_tbi = work_dir + "/out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 2
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


# Remark: since the input comes from $IntraSampleMerge$, which overwrites GTs,
# every call in the output of this procedure has at least one sample with
# GT=0/1, and every sample has GT \in {0/1, ./.}.
#
# Performance on 1074 AoU 8x samples:
# COMMAND           RUNTIME     N_CPUS      MAX_RSS
# bcftools merge    2.5h        2           140G
# bcftools norm     40m         3           80G
# awk               2h          0.5         75G
# bgzip             10m         2           23M
#
task InterSampleMerge {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
        Int compression_level = 1
        Int ram_gb = 200
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(input_vcf_gz, "GB")) + 100
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_gb} - 2 ))
        
        # Merging all intra-sample VCFs
        INPUT_FILES=~{sep=',' input_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        for INPUT_FILE in ${INPUT_FILES}; do
            echo ${INPUT_FILE} >> list.txt
        done
        ls -laht; tree; df -h
        # $--info-rules -$ disables default rules, and it is used just to avoid
        # the following error:
        # Only fixed-length vectors are supported with -i sum:AC
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --info-rules - --file-list list.txt --output-type z > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        ls -laht; tree; df -h
        
        # Removing multiallelic records, if any are generated during the merge.
        # This is just an extra safeguard and might be dropped if $bcftools
        # merge --merge none$ always behaves correctly.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z tmp1.vcf.gz > tmp2.vcf.gz
        tabix -f tmp2.vcf.gz
        rm -f tmp1.vcf.gz*
        ls -laht; tree; df -h
        
        # Restoring symbolic ALTs to their original states
        bcftools view --header-only tmp2.vcf.gz > merged.vcf
        ${TIME_COMMAND} bcftools view --no-header tmp2.vcf.gz | awk '{ \
            tag="artificial"; \
            if (substr($5,6,length(tag))==tag) $5=substr($5,1,4) ">"; \
            printf("%s",$1); \
            for (i=2; i<=NF; i++) printf("\t%s",$i); \
            printf("\n"); \
        }' >> merged.vcf
        rm -f tmp2.vcf.gz*
        ${TIME_COMMAND} bgzip --threads ${N_THREADS} --compress-level ~{compression_level} merged.vcf
        tabix -f merged.vcf.gz
        ls -laht; tree; df -h
    >>>
    
    output {
        File output_vcf_gz = work_dir + "/merged.vcf.gz"
        File output_tbi = work_dir + "/merged.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 4
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
