version 1.0
import "hapestry_annotate.wdl" as hapestry_annotate


struct RuntimeAttributes {
    Int? cpu
    Int? command_mem_gb
    Int? additional_mem_gb
    Int? disk_size_gb
    Int? boot_disk_size_gb
    Boolean? use_ssd
    Int? preemptible
    Int? max_retries
}


# Define the task
task chunk_vcf {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File? confident_bed

        # Hapestry specific args
        Int interval_max_length = 50000
        Int flank_length = 200
        Int min_sv_length = 20
        Int n_chunks = 32
        File tandems_bed

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    String docker_dir = "/hapestry"
    String output_dir = "output"

    command <<<
        set -eoxu pipefail

        mkdir ~{output_dir}

        # get the directory of the vcf
        vcf_dir=$(dirname ~{vcf_gz})

        # get the tbi directory
        tbi_dir=$(dirname ~{vcf_gz_tbi})

        # if the tbi_dr and the vcf_dir are not the same then move the tbi into the vcf directory
        if [ "$vcf_dir" != "$tbi_dir" ]; then
            mv ~{vcf_gz_tbi} ${vcf_dir}
        fi

        if [ -s ~{monitoring_script} ]; then
          bash ~{monitoring_script} > monitoring.log &
        fi

        if ~{defined(confident_bed)}
        then
            # use bcftools to subset the vcf by the confident bed
            bcftools view -T ~{confident_bed} ~{vcf_gz} -Ov -o confident.vcf

            # convert to bgzipped vcf and overwrite the input VCF
            bcftools view -Oz -o ~{vcf_gz} confident.vcf

            # index the vcf
            bcftools index -t ~{vcf_gz}
        fi

        # create a variable for the non-gz vcf
        vcf=$(basename ~{vcf_gz} .gz)

        # unzip the VCF for hapestry
        bcftools view -Ov -o ${vcf} ~{vcf_gz}

        ~{docker_dir}/sv_merge/build/find_windows \
        --output_dir ~{output_dir}/run/ \
        --n_chunks ~{n_chunks} \
        --vcf ${vcf} \
        --tandems ~{tandems_bed} \
        --interval_max_length ~{interval_max_length} \
        --min_sv_length ~{min_sv_length} \
        --flank_length ~{flank_length}

        tree ~{output_dir}

        # use each of the generated BEDs to subset the VCF
        for file in ~{output_dir}/run/*; do
            [ -e "$file" ] || continue

            # if the file is not of the format windows_[numeric].bed, using regex to identify the pattern, skip it
            if [[ ! $(basename ${file}) =~ ^windows_[0-9]+\.bed$ ]]; then
                continue
            fi

            echo "processing ${file}"
            bcftools view -T ${file} -Oz -o "~{output_dir}/$(basename ${file}).vcf.gz" ~{vcf_gz}
            bcftools index -t -o "~{output_dir}/$(basename ${file}).vcf.gz.tbi" "~{output_dir}/$(basename ${file}).vcf.gz"
        done

        tree ~{output_dir}
        >>>

    parameter_meta {
        interval_max_length: "Maximum length of each window evaluated"
        flank_length: "Length of flanking sequence to include in each window"
        tandems_bed: "BED file of tandem repeats"
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        Array[File] chunked_vcfs = glob("output/*.vcf.gz")
        Array[File] chunked_tbis = glob("output/*.vcf.gz.tbi")
    }
}


# A task that takes Array[File] VCF and concatenates them
task bcftools_concat{
    input {
        Array[File] vcf_gz
        Array[File] vcf_gz_tbi

        String docker = "staphb/bcftools:1.20"

        RuntimeAttributes runtime_attributes = {}
    }

    command <<<
        set -eou pipefail

        # Concatenate the VCFs
        bcftools concat -a -D -Oz -o concatenated.vcf.gz ~{sep=" " vcf_gz}
        bcftools sort -Oz -o concatenated_sorted.vcf.gz concatenated.vcf.gz
        bcftools index -t concatenated_sorted.vcf.gz
    >>>

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File concatenated_vcf_gz = "concatenated_sorted.vcf.gz"
        File concatenated_vcf_gz_tbi = "concatenated_sorted.vcf.gz.tbi"
    }
}


workflow hapestry_annotate_scattered {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File confident_bed

        # Hapestry specific args
        Int interval_max_length = 50000
        Int flank_length = 200
        Int n_threads
        Int n_chunks
        File tandems_bed
        File reference_fa
        File haps_vs_ref_csv
        Boolean force_unique_reads = false
        Boolean bam_not_hardclipped = false
        String annotation_label = "HAPESTRY_REF"

        String docker
        File? monitoring_script

        RuntimeAttributes? annotate_runtime_attributes
    }

    parameter_meta {
        interval_max_length: "Maximum length of each window evaluated"
        flank_length: "Length of flanking sequence to include in each window"
        min_sv_length: "Minimum SV length to consider"
        n_chunks: "Number of chunks to split the VCF into, and subsequently the number of workers"
        n_threads: "Maximum number of threads to use"
        tandems_bed: "BED file of tandem repeats"
        reference_fa: "Reference fasta file"
        haps_vs_ref_csv: "CSV file of haplotype vs reference BAMs"
        force_unique_reads: "Force unique aligned sequence names among multiple BAMs to prevent collisions"
        bam_not_hardclipped: "If the bam is GUARANTEED not to contain any hardclips, use this flag to trigger much simpler/faster fetching process"
        annotation_label: "Name to give the INFO field in the VCF for annotations, usually upper case"
    }

    call chunk_vcf{
        input:
            vcf_gz = vcf_gz,
            vcf_gz_tbi = vcf_gz_tbi,
            confident_bed = confident_bed,
            interval_max_length = interval_max_length,
            flank_length = flank_length,
            min_sv_length = min_sv_length,
            tandems_bed = tandems_bed,
            n_chunks = n_chunks,
            docker = docker,
            monitoring_script = monitoring_script
    }

    Array[Pair[File,File]] items = zip(chunk_vcf.chunked_vcfs, chunk_vcf.chunked_tbis)

    scatter (x in items){
        call hapestry_annotate.annotate as scattered_annotate {
            input:
                vcf_gz = x.left,
                vcf_gz_tbi = x.right,
                bam_not_hardclipped = bam_not_hardclipped,
                interval_max_length = interval_max_length,
                flank_length = flank_length,
                n_threads = n_threads,
                tandems_bed = tandems_bed,
                reference_fa = reference_fa,
                haps_vs_ref_csv = haps_vs_ref_csv,
                force_unique_reads = force_unique_reads,
                annotation_label = annotation_label,
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = annotate_runtime_attributes
        }
    }

    # Call concat
    call bcftools_concat {
        input:
            vcf_gz = scattered_annotate.annotated_vcf_gz,
            vcf_gz_tbi = scattered_annotate.annotated_vcf_gz_tbi,
            docker = docker
    }

    output {
        File annotated_vcf_gz = bcftools_concat.concatenated_vcf_gz
        File annotated_vcf_gz_tbi = bcftools_concat.concatenated_vcf_gz_tbi
    }
}
