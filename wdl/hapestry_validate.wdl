version 1.0

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
task validate {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File training_resource_bed

        # Hapestry specific args
        Float? min_score = 0.9
        Int? interval_max_length = 50000
        Int? flank_length = 200
        Int n_threads
        File tandems_bed
        File reference_fa
        File haps_vs_ref_csv
        Boolean force_unique_reads = false
        String? annotation_label = "HAPESTRY_REF"

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    String docker_dir = "/hapestry"
    String output_dir = "output/"

    command <<<
        set -eou pipefail

        if [ -s ~{monitoring_script} ]; then
          bash ~{monitoring_script} > monitoring.log &
        fi

        # use bcftools to subset the vcf by the training resource bed
        bcftools view -T ~{training_resource_bed} ~{vcf_gz} -Ov -o training.vcf

        ~{docker_dir}/sv_merge/build/annotate \
        --output_dir ~{output_dir}/run/ \
        --bam_csv ~{haps_vs_ref_csv} \
        --vcf training.vcf \
        --tandems ~{tandems_bed} \
        --ref ~{reference_fa} \
        --interval_max_length ~{interval_max_length} \
        --flank_length ~{flank_length} \
        --n_threads ~{n_threads} \
        --label ~{annotation_label} ~{if force_unique_reads then "--force_unique_reads" else ""}

        # use bcftools view -i to filter the true positive sites according to the annotation_label + _MAX > min_score
        bcftools view -i "INFO/~{annotation_label}_MAX > ~{min_score}" ~{output_dir}/run/annotated.vcf -Oz -o ~{output_dir}/run/annotated_tp.vcf.gz

        # index the training vcf
        bcftools index -t ~{output_dir}/run/annotated_tp.vcf.gz
    >>>

    parameter_meta {
        min_score: "Minimum identity observed spanning a variant to consider a true positive w.r.t. assembly alignments"
        interval_max_length: "Maximum length of each window evaluated"
        flank_length: "Length of flanking sequence to include in each window"
        n_threads: "Maximum number of threads to use"
        tandems_bed: "BED file of tandem repeats"
        reference_fa: "Reference fasta file"
        haps_vs_ref_csv: "CSV file of haplotype vs reference BAMs"
        force_unique_reads: "Force unique aligned sequence names among multiple BAMs to prevent collisions"
        annotation_label: "Name to give the INFO field in the VCF for annotations, usually upper case"
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
        File training_vcf_gz = output_dir + "/run/annotated_tp.vcf.gz"
        File training_vcf_gz_tbi = output_dir + "/run/annotated_tp.vcf.gz.tbi"
        File? monitoring_log = "monitoring.log"
    }
}


# Create main workflow which calls the IdentifyTrainingSites task
workflow hapestry_validate {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File training_resource_bed

        # Hapestry specific args
        Float? min_score = 0.9
        Int? interval_max_length = 50000
        Int? flank_length = 200
        Int n_threads
        File tandems_bed
        File reference_fa
        File haps_vs_ref_csv
        Boolean force_unique_reads = false
        String? annotation_label = "HAPESTRY_REF"

        String docker
        File? monitoring_script

        RuntimeAttributes? validate_runtime_attributes
    }

    parameter_meta {
        min_score: "Minimum identity observed spanning a variant to consider a true positive w.r.t. assembly alignments"
        interval_max_length: "Maximum length of each window evaluated"
        flank_length: "Length of flanking sequence to include in each window"
        n_threads: "Maximum number of threads to use"
        tandems_bed: "BED file of tandem repeats"
        reference_fa: "Reference fasta file"
        haps_vs_ref_csv: "CSV file of haplotype vs reference BAMs"
        force_unique_reads: "Force unique aligned sequence names among multiple BAMs to prevent collisions"
        annotation_label: "Name to give the INFO field in the VCF for annotations, usually upper case"
    }

    call validate {
        input:
            vcf_gz = vcf_gz,
            vcf_gz_tbi = vcf_gz_tbi,
            training_resource_bed = training_resource_bed,
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
            runtime_attributes = validate_runtime_attributes
    }

    output {
        File training_vcf_gz = validate.training_vcf_gz
        File training_vcf_gz_tbi = validate.training_vcf_gz_tbi
        File? monitoring_log = validate.monitoring_log
    }
}


