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
task annotate {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File? confident_bed

        # Hapestry specific args
        Int? interval_max_length = 50000
        Int? flank_length = 200
        Int n_threads
        File tandems_bed
        File reference_fa
        File haps_vs_ref_csv
        Boolean force_unique_reads = false
        Boolean bam_not_hardclipped = false
        String? annotation_label = "HAPESTRY_REF"

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    String docker_dir = "/hapestry"
    String output_dir = "output/"

    command <<<
        git --no-pager --git-dir ~{docker_dir}/sv_merge/.git log --decorate=short --pretty=oneline | head -n 1

        set -eou pipefail

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

        # use bcftools to subset the vcf by the confident bed, only if the bed is defined
        if ~{defined(confident_bed)}; then
            bcftools view -R ~{confident_bed} ~{vcf_gz} -Ov -o confident.vcf
        else
            bcftools view -Ov ~{vcf_gz} -o confident.vcf
        fi

        ~{docker_dir}/sv_merge/build/annotate \
        --output_dir ~{output_dir}/run/ \
        --bam_csv ~{haps_vs_ref_csv} \
        --vcf confident.vcf \
        --tandems ~{tandems_bed} \
        --ref ~{reference_fa} \
        --interval_max_length ~{interval_max_length} \
        --flank_length ~{flank_length} \
        --n_threads ~{n_threads} \
        --label ~{annotation_label} ~{if force_unique_reads then "--force_unique_reads" else ""} ~{if bam_not_hardclipped then "--bam_not_hardclipped" else ""}

        # convert to bgzipped vcf
        bcftools view -Oz -o ~{output_dir}/run/annotated.vcf.gz ~{output_dir}/run/annotated.vcf

        # index the vcf
        bcftools index -t ~{output_dir}/run/annotated.vcf.gz

        # tarball the output directory -- WARNING omitted because tar takes 2 hours on these outputs...
        # tar -cvzf ~{output_dir}/hapestry.tar.gz ~{output_dir}/run
    >>>

    parameter_meta {
        interval_max_length: "Maximum length of each window evaluated"
        flank_length: "Length of flanking sequence to include in each window"
        n_threads: "Maximum number of threads to use"
        tandems_bed: "BED file of tandem repeats"
        confident_bed: "BED file of regions to be included, if not provided, all variants will be annotated"
        reference_fa: "Reference fasta file"
        haps_vs_ref_csv: "CSV file of haplotype vs reference BAMs"
        force_unique_reads: "Force unique aligned sequence names among multiple BAMs to prevent collisions"
        bam_not_hardclipped: "If the bam is GUARANTEED not to contain any hardclips, use this flag to trigger much simpler/faster fetching process"
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
        File annotated_vcf_gz = output_dir + "/run/annotated.vcf.gz"
        File annotated_vcf_gz_tbi = output_dir + "/run/annotated.vcf.gz.tbi"
#        File hapestry_annotate_data = output_dir + "/hapestry.tar.gz"
        File? monitoring_log = "monitoring.log"
    }
}


workflow hapestry_annotate {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File confident_bed

        # Hapestry specific args
        Int? interval_max_length = 50000
        Int? flank_length = 200
        Int n_threads
        File tandems_bed
        File reference_fa
        File haps_vs_ref_csv
        Boolean force_unique_reads = false
        Boolean bam_not_hardclipped = false
        String? annotation_label = "HAPESTRY_REF"

        String docker
        File? monitoring_script

        RuntimeAttributes? annotate_runtime_attributes
    }

    parameter_meta {
        interval_max_length: "Maximum length of each window evaluated"
        flank_length: "Length of flanking sequence to include in each window"
        n_threads: "Maximum number of threads to use"
        tandems_bed: "BED file of tandem repeats"
        reference_fa: "Reference fasta file"
        haps_vs_ref_csv: "CSV file of haplotype vs reference BAMs"
        force_unique_reads: "Force unique aligned sequence names among multiple BAMs to prevent collisions"
        bam_not_hardclipped: "If the bam is GUARANTEED not to contain any hardclips, use this flag to trigger much simpler/faster fetching process"
        annotation_label: "Name to give the INFO field in the VCF for annotations, usually upper case"
    }

    call annotate {
        input:
            vcf_gz = vcf_gz,
            vcf_gz_tbi = vcf_gz_tbi,
            confident_bed = confident_bed,
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

    output {
        File annotated_vcf_gz = annotate.annotated_vcf_gz
        File annotated_vcf_gz_tbi = annotate.annotated_vcf_gz_tbi
#        File hapestry_annotate_data = annotate.hapestry_annotate_data
        File? monitoring_log = annotate.monitoring_log
    }
}
