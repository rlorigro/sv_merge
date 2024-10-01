version 1.0
import "hapestry_merge.wdl" as hapestry_merge


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
        File? windows_override_bed
        File? reference_fa

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

        # Kind of circular to provide windows to "find_windows" app, but if we dont then we have to reproduce a lot of
        # clean up and chunking separately in bash, which is a waste of time.
        windows_arg=""
        if ~{defined(windows_override_bed)}
        then
            windows_arg="--windows ~{windows_override_bed}"
        fi

        reference_arg=""
        if ~{defined(reference_fa)}
        then
            reference_arg="--ref_fasta ~{reference_fa}"
        fi

        ~{docker_dir}/sv_merge/build/find_windows \
        --output_dir ~{output_dir}/run/ \
        --n_chunks ~{n_chunks} \
        --vcf ${vcf} \
        --tandems ~{tandems_bed} \
        --interval_max_length ~{interval_max_length} \
        --min_sv_length ~{min_sv_length} \
        --flank_length ~{flank_length} ${windows_arg} ${reference_arg}

        tree ~{output_dir}

        # use each of the generated BEDs to subset the VCF
        for file in ~{output_dir}/run/*; do
            [ -e "$file" ] || continue

            # if the file is not of the format windows_[numeric]_unflanked.bed, using regex to identify the pattern, skip it
            if [[ ! $(basename ${file}) =~ ^windows_[0-9]+_unflanked\.bed$ ]]; then
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
        windows_override_bed: "BED file of windows (NOT including flank regions!) to be used for merging. Flanks will be added to the windows provided."
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


# Task to combine the chunks
task concat_vcfs{
    input {
        Array[File] sequence_tarballs

        String docker = "staphb/bcftools:1.20"

        RuntimeAttributes runtime_attributes = {}
    }

    command <<<
    set -eoxu pipefail

    # Temporary directory for extracted VCF files
    temp_dir=$(mktemp -d)
    echo "Temporary directory created at: $temp_dir"

    # Array to hold paths of extracted VCF files
    vcf_files=()

    i=0

    # Loop through each tar.gz file provided as an argument
    for archive in ~{sep=' ' sequence_tarballs}; do
        echo "Processing archive: $archive"

        # find the path of the VCF file in the tarball and don't give error exit code if not found
        vcf_path=$(tar -tzf "$archive" | grep -oP '.*\merged.vcf$') || true

        echo $vcf_path

        # check if 'vcf' in the path
        if [[ $vcf_path != *"vcf"* ]]; then
            echo "No VCF file found in archive: $archive"
            continue
        fi

        # Extract the VCF file
        tar -xzf "$archive" -C "$temp_dir" "$vcf_path"

        # Pick a new name for the extracted VCF file and then rename it to avoid conflicts
        new_vcf_path="$temp_dir/extracted_vcf_${i}.vcf.gz"

        echo "new path: $new_vcf_path"

        bcftools view -Oz -o "$new_vcf_path" "$temp_dir/$vcf_path"
        bcftools index -t "$new_vcf_path"

        # Add the path of the extracted VCF to the array
        vcf_files+=("$new_vcf_path")
        i=$((i+1))
    done

    bcftools concat -a -D -Oz -o concatenated.vcf.gz "${vcf_files[@]}"
    bcftools sort -Oz -o concatenated_sorted.vcf.gz concatenated.vcf.gz
    bcftools index -t concatenated_sorted.vcf.gz

    >>>

    output {
        File concatenated_vcf = "concatenated_sorted.vcf.gz"
        File concatenated_vcf_tbi = "concatenated_sorted.vcf.gz.tbi"
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
}


# Task to combine the chunks
task concat_beds{
    input {
        Array[File] bed_tarballs

        String docker = "staphb/bedtools:2.31.1"

        RuntimeAttributes runtime_attributes = {}
    }

    command <<<
    set -eoxu pipefail

    # iterate the tarballs. each one will contain BED files like so:
    # - windows.bed
    # - windows_unflanked.bed
    # - windows_failed.bed
    # - windows_omitted.bed
    # for each tarball, extract the BEDs and then append them to a growing BED for each filename

    # Temporary directory for extracted BED files
    temp_dir=$(mktemp -d)
    echo "Temporary directory created at: $temp_dir"

    # Iterate the tarballs and concatenate
    for archive in ~{sep=' ' bed_tarballs}; do
        echo "Processing archive: $archive"

        # extract the BED files
        tar -xzf "$archive" -C "$temp_dir"

        # append the BEDs to the growing BEDs
        for file in $(find "$temp_dir" -type f); do
            [ -e "$file" ] || continue

            echo "processing ${file}"
            # Check if the item is a file before processing
            if [ -f "$file" ]; then
                echo "is file: ${file}"
                cat ${file} >> $(basename ${file})
            fi
        done

        # remove the extracted files
        rm -rf "$temp_dir"/*
    done

    # tarball the BEDs
    tar -cvzf beds.tar.gz ./*.bed

    >>>

    output {
        File beds_tarball = "beds.tar.gz"
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
}


workflow hapestry_merge_scattered {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File confident_bed
        File? windows_override_bed

        Int interval_max_length = 50000
        Int flank_length = 200
        Int min_sv_length = 20
        Int graphaligner_timeout = 120
        Int solver_timeout = 30*60
        Float min_read_hap_identity = 0.5
        Float d_weight = 1.0
        Int n_threads
        Int n_chunks
        File tandems_bed
        File reference_fa
        File haps_vs_ref_csv
        Boolean force_unique_reads = false
        Boolean bam_not_hardclipped = false
        Boolean skip_solve = false
        Boolean quadratic_objective = false
        Boolean rescale_weights = false

        String docker
        File? monitoring_script

        RuntimeAttributes? merge_runtime_attributes
    }

    parameter_meta {
        windows_override_bed: "BED file of windows (NOT including flank regions!) to be used for merging. Flanks will be added to the windows provided."
        bam_not_hardclipped: "If the bam is GUARANTEED not to contain any hardclips, use this flag to trigger much simpler/faster fetching process"
        flank_length: "Length of flanking sequence to include in each window"
        force_unique_reads: "Force unique aligned sequence names among multiple BAMs to prevent collisions"
        graphaligner_timeout: "Timeout for graphaligner in seconds"
        solver_timeout: "Timeout for solver in seconds"
        haps_vs_ref_csv: "CSV file of haplotype vs reference BAMs"
        interval_max_length: "Maximum length of each window evaluated"
        min_read_hap_identity: "Minimum identity between read and haplotype to consider as input to optimizer"
        d_weight: "Scaling factor for the distance term in the objective function, larger than 1.0 gives greater priority to minimizing edit distance"
        min_sv_length: "Only variants that affect at least this number of bps are merged. Shorter variants are used to build graphs and haplotypes, but they are not merged or printed in output."
        n_chunks: "Number of chunks to split the VCF into, and subsequently the number of workers"
        n_threads: "Maximum number of threads to use"
        reference_fa: "Reference fasta file"
        skip_solve: "Skip the solve step, only generate input CSV for the solve step"
        quadratic_objective: "Use quadratic objective which finds the normalized square distance from the utopia point"
        rescale_weights: "Use quadratic difference-from-best scaling for weights"
        tandems_bed: "BED file of tandem repeats"
    }

    call chunk_vcf {
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
            monitoring_script = monitoring_script,
            reference_fa = reference_fa,
            windows_override_bed = windows_override_bed
    }

    Array[Pair[File,File]] items = zip(chunk_vcf.chunked_vcfs, chunk_vcf.chunked_tbis)

    scatter (x in items) {
        call hapestry_merge.merge as scattered_merge {
            input:
                vcf_gz = x.left,
                vcf_gz_tbi = x.right,
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
                quadratic_objective = quadratic_objective,
                rescale_weights = rescale_weights,
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = merge_runtime_attributes,
                confident_bed = confident_bed
        }
    }

    call concat_vcfs {
        input:
            sequence_tarballs = scattered_merge.sequence_data_tarball
    }

    call concat_beds {
        input:
            bed_tarballs = scattered_merge.beds_tarball
    }

    output {
        Array[File] non_sequence_data_tarball = scattered_merge.non_sequence_data_tarball
        Array[File] sequence_data_tarball = scattered_merge.sequence_data_tarball
        File hapestry_vcf = concat_vcfs.concatenated_vcf
        File hapestry_vcf_tbi = concat_vcfs.concatenated_vcf_tbi
        File hapestry_beds_tarball = concat_beds.beds_tarball
    }
}
