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
task merge {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File? confident_bed

        # Hapestry specific args
        Int interval_max_length = 50000
        Int min_sv_length = 20
        Int flank_length = 200
        Int graphaligner_timeout = 120
        Int solver_timeout = 30*60
        Float min_read_hap_identity = 0.5
        Float d_weight = 1.0
        Int n_threads
        File tandems_bed
        File reference_fa
        File haps_vs_ref_csv
        Boolean force_unique_reads = false
        Boolean bam_not_hardclipped = false
        Boolean skip_solve = false
        Boolean rescale_weights = false

        String docker = "fcunial/hapestry:merge"
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    String docker_dir = "/hapestry"
    String output_dir = "output/"

    command <<<
        git --no-pager --git-dir ~{docker_dir}/sv_merge/.git log --decorate=short --pretty=oneline | head -n 1

        set -eoxu pipefail

        # get the directory of the vcf
        vcf_dir=$(dirname ~{vcf_gz})

        # get the tbi directory
        tbi_dir=$(dirname ~{vcf_gz_tbi})

        # if the tbi_dir and the vcf_dir are not the same then move the tbi into the vcf directory
        if [ "$vcf_dir" != "$tbi_dir" ]; then
            mv ~{vcf_gz_tbi} ${vcf_dir}
        fi

        if [ -s ~{monitoring_script} ]; then
          bash ~{monitoring_script} > monitoring.log &
        fi

        # use bcftools to subset the vcf by the confident bed, only if the bed is defined
        if ~{defined(confident_bed)}; then
            bcftools view -T ~{confident_bed} ~{vcf_gz} -Ov -o confident.vcf
        else
            bcftools view -Ov ~{vcf_gz} -o confident.vcf
        fi

        ~{docker_dir}/sv_merge/build/hapestry \
        --output_dir ~{output_dir}/run/ \
        --bam_csv ~{haps_vs_ref_csv} \
        --vcf confident.vcf \
        --tandems ~{tandems_bed} \
        --ref ~{reference_fa} \
        --interval_max_length ~{interval_max_length} \
        --min_sv_length ~{min_sv_length} \
        --flank_length ~{flank_length} \
        --graphaligner_timeout ~{graphaligner_timeout} \
        --solver_timeout ~{solver_timeout} \
        --min_read_hap_identity ~{min_read_hap_identity} \
        --d_weight ~{d_weight} \
        --n_threads ~{n_threads} \
        ~{if force_unique_reads then "--force_unique_reads" else ""} \
        ~{if bam_not_hardclipped then "--bam_not_hardclipped" else ""} \
        ~{if skip_solve then "--skip_solve" else ""} \
        ~{if rescale_weights then "--rescale_weights" else ""}

        # tarball only the csv files in the output subdirectories
        find ~{output_dir}/run/ \( -name "*.csv" -o -name "*.txt" \) -exec tar -cvzf ~{output_dir}/non_sequence_data.tar.gz {} +
        find ~{output_dir}/run/ \( -name "*.fasta" -o -name "*.gfa" -o -name "*.gaf" -o -name "*.vcf" \) -exec tar -cvzf ~{output_dir}/sequence_data.tar.gz {} +

        # if the outputs are empty, create empty placeholders
        if [ ! -s ~{output_dir}/non_sequence_data.tar.gz ]; then
            tar -cvzf ~{output_dir}/non_sequence_data.tar.gz --files-from /dev/null
        fi
        if [ ! -s ~{output_dir}/sequence_data.tar.gz ]; then
            tar -cvzf ~{output_dir}/sequence_data.tar.gz --files-from /dev/null
        fi

        # tarball just the BED files in the top level output directory
        tar -cvzf ~{output_dir}/beds.tar.gz ~{output_dir}/run/*.bed
    >>>

    parameter_meta {
        bam_not_hardclipped: "If the bam is GUARANTEED not to contain any hardclips, use this flag to trigger much simpler/faster fetching process"
        confident_bed: "BED file of regions to be included, if not provided, all variants will be merged"
        flank_length: "Length of flanking sequence to include in each window"
        force_unique_reads: "Force unique aligned sequence names among multiple BAMs to prevent collisions"
        graphaligner_timeout: "Timeout for graphaligner in seconds"
        solver_timeout: "Timeout for solver in seconds"
        haps_vs_ref_csv: "CSV file of haplotype vs reference BAMs"
        interval_max_length: "Maximum length of each window evaluated"
        min_read_hap_identity: "Minimum identity between read and haplotype to consider as input to optimizer"
        d_weight: "Scaling factor for the D term in the optimizer, greater than 1.0 will prioritize minimizing edit distance"
        min_sv_length: "Minimum length of SVs to consider"
        n_threads: "Maximum number of threads to use"
        reference_fa: "Reference fasta file"
        skip_solve: "Skip the solve step, only generate input CSV for the solve step"
        rescale_weights: "Use quadratic difference-from-best scaling for weights"
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
        File non_sequence_data_tarball = output_dir + "/non_sequence_data.tar.gz"
        File sequence_data_tarball = output_dir + "/sequence_data.tar.gz"
        File beds_tarball = output_dir + "/beds.tar.gz"
        File? monitoring_log = "monitoring.log"
    }
}


workflow hapestry_merge {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File confident_bed

        # Hapestry specific args
        Int interval_max_length = 50000
        Int flank_length = 200
        Int min_sv_length = 20
        Int graphaligner_timeout = 120
        Int solver_timeout = 30*60
        Float min_read_hap_identity = 0.5
        Float d_weight = 1.0
        Int n_threads
        File tandems_bed
        File reference_fa
        File haps_vs_ref_csv
        Boolean force_unique_reads = false
        Boolean bam_not_hardclipped = false
        Boolean skip_solve = false
        Boolean rescale_weights = false

        String docker
        File? monitoring_script

        RuntimeAttributes? merge_runtime_attributes
    }

    parameter_meta {
        bam_not_hardclipped: "If the bam is GUARANTEED not to contain any hardclips, use this flag to trigger much simpler/faster fetching process"
        confident_bed: "BED file of regions to be included, if not provided, all variants will be merged"
        flank_length: "Length of flanking sequence to include in each window"
        force_unique_reads: "Force unique aligned sequence names among multiple BAMs to prevent collisions"
        graphaligner_timeout: "Timeout for graphaligner in seconds"
        solver_timeout: "Timeout for solver in seconds"
        haps_vs_ref_csv: "CSV file of haplotype vs reference BAMs"
        interval_max_length: "Maximum length of each window evaluated"
        min_read_hap_identity: "Minimum identity between read and haplotype to consider as input to optimizer"
        d_weight: "Scaling factor for the D term in the optimizer, greater than 1.0 will prioritize minimizing edit distance"
        min_sv_length: "Minimum length of SVs to consider"
        n_threads: "Maximum number of threads to use"
        reference_fa: "Reference fasta file"
        skip_solve: "Skip the solve step, only generate input CSV for the solve step"
        rescale_weights: "Use quadratic difference-from-best scaling for weights"
        tandems_bed: "BED file of tandem repeats"
    }

    call merge {
        input:
            vcf_gz = vcf_gz,
            vcf_gz_tbi = vcf_gz_tbi,
            confident_bed = confident_bed,
            bam_not_hardclipped = bam_not_hardclipped,
            interval_max_length = interval_max_length,
            flank_length = flank_length,
            min_sv_length = min_sv_length,
            graphaligner_timeout = graphaligner_timeout,
            solver_timeout = solver_timeout,
            min_read_hap_identity = min_read_hap_identity,
            d_weight = d_weight,
            n_threads = n_threads,
            tandems_bed = tandems_bed,
            reference_fa = reference_fa,
            haps_vs_ref_csv = haps_vs_ref_csv,
            force_unique_reads = force_unique_reads,
            skip_solve = skip_solve,
            rescale_weights = rescale_weights,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = merge_runtime_attributes
    }

    output {
        File non_sequence_data_tarball = merge.non_sequence_data_tarball
        File sequence_data_tarball = merge.sequence_data_tarball
        File beds_tarball = merge.beds_tarball
        File? monitoring_log = merge.monitoring_log
    }
}
