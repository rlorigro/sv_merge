version 1.0


# 
#
workflow GraphEvaluation {
    input {
        Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
        Int interval_max_length
        Int flank_length
        File reference_fa
        File tandems_bed
        Array[File] evaluation_bed
        Array[String] intersample_vcf_gz
        File haps_vs_chm13_csv
    }
    parameter_meta {
        intersample_vcf_gz: "List of remote VCF.GZ addresses. The evaluation builds haplotype clusters using the first file."
        haps_vs_chm13_csv: "List of remote BAM addresses"
    }
    
    scatter(chromosome in chromosomes) {
        call EvaluateChromosome {
            input:
                chromosome = chromosome,
                interval_max_length = interval_max_length,
                flank_length = flank_length,
                reference_fa = reference_fa,
                tandems_bed = tandems_bed,
                evaluation_bed = evaluation_bed,
                intersample_vcf_gz = intersample_vcf_gz,
                haps_vs_chm13_csv = haps_vs_chm13_csv
        }
    }
    
    output {
        Array[File] analysis = EvaluateChromosome.analysis
        Array[File] evaluation = EvaluateChromosome.evaluation
    }
}


# Evaluates multiple VCFs on a given chromosome, stratifying the analysis over
# several BED files.
#
task EvaluateChromosome {
    input {
        String chromosome
        Int interval_max_length
        Int flank_length
        File reference_fa
        File tandems_bed
        Array[File] evaluation_bed
        Array[String] intersample_vcf_gz
        File haps_vs_chm13_csv
    }
    parameter_meta {
    }
    
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
#        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        # TEMPORARILY FORCE 1 THREAD
        N_THREADS=1

        # Downloading all the calls in $chromosome$.
        export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
        LIST_FILE=~{write_lines(intersample_vcf_gz)}
        VCFS=""; CLUSTER_BY=""; TOOLS=""; LINE="0";
        while read ADDRESS; do
            TOOL_ID=$(basename ${ADDRESS} .vcf.gz)
            bcftools view --output-type z ${ADDRESS} ~{chromosome} > tmp.vcf.gz
            tabix tmp.vcf.gz
            bcftools norm --multiallelics - --output-type v tmp.vcf.gz > ${TOOL_ID}.vcf
            rm -f tmp.vcf.gz
            LINE=$(( ${LINE} + 1 ))
            if [[ ${LINE} -eq 1 ]]; then
                VCFS=${TOOL_ID}.vcf
                CLUSTER_BY=${TOOL_ID}.vcf
                TOOLS=${TOOL_ID}
            else
                VCFS="${VCFS},${TOOL_ID}.vcf"
                TOOLS="${TOOLS} ${TOOL_ID}"
            fi
        done < ${LIST_FILE}
        
        # Evaluating
        EVALUATION_BEDS=~{sep=',' evaluation_bed}
        EVALUATION_BEDS=$(echo ${EVALUATION_BEDS} | tr ',' ' ')
        EVALUATION_NAME="~{chromosome}_evaluation"
        ANALYSIS_NAME="~{chromosome}_analysis"
        ${TIME_COMMAND} ~{docker_dir}/sv_merge/build/evaluate \
            --n_threads ${N_THREADS} \
            --output_dir ./${EVALUATION_NAME} \
            --bam_csv ~{haps_vs_chm13_csv} \
            --vcfs ${VCFS} \
            --cluster_by ${CLUSTER_BY} \
            --tandems ~{tandems_bed} \
            --ref ~{reference_fa} \
            --interval_max_length ~{interval_max_length} \
            --flank_length ~{flank_length} \
            --debug
        mkdir ./${ANALYSIS_NAME}
        ${TIME_COMMAND} ~{docker_dir}/sv_merge/build/analyze_evaluation \
            --input_dir ./${EVALUATION_NAME} \
            --output_dir ./${ANALYSIS_NAME} \
            --tools ${TOOLS} \
            --beds ${EVALUATION_BEDS}
        
        export GZIP=-1
        ${TIME_COMMAND} tar -czf ${ANALYSIS_NAME}.tar.gz ./${ANALYSIS_NAME}
        ${TIME_COMMAND} tar -czf ${EVALUATION_NAME}.tar.gz --exclude='*.fasta' --exclude='*.fa' ./${EVALUATION_NAME}
    >>>
    
    output {
        File analysis = work_dir + "/" + chromosome + "_analysis.tar.gz"
        File evaluation = work_dir + "/" + chromosome + "_evaluation.tar.gz"
    }
    runtime {
        docker: "fcunial/hapestry"
        cpu: 32
        memory: 128
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}
