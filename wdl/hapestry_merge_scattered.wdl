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


workflow hapestry_merge_scattered {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File confident_bed

        # Hapestry specific args
        Int interval_max_length = 50000
        Int flank_length = 200
        Int min_sv_length = 20
        Int graphaligner_timeout = 120
        Float min_read_hap_identity = 0.5
        Int n_threads
        Int n_chunks
        File tandems_bed
        File reference_fa
        File haps_vs_ref_csv
        Boolean force_unique_reads = false
        Boolean bam_not_hardclipped = false
        Boolean skip_solve = false

        String docker
        File? monitoring_script

        RuntimeAttributes? merge_runtime_attributes
    }

    parameter_meta {
        bam_not_hardclipped: "If the bam is GUARANTEED not to contain any hardclips, use this flag to trigger much simpler/faster fetching process"
        flank_length: "Length of flanking sequence to include in each window"
        force_unique_reads: "Force unique aligned sequence names among multiple BAMs to prevent collisions"
        graphaligner_timeout: "Timeout for graphaligner in seconds"
        haps_vs_ref_csv: "CSV file of haplotype vs reference BAMs"
        interval_max_length: "Maximum length of each window evaluated"
        min_read_hap_identity: "Minimum identity between read and haplotype to consider as input to optimizer"
        min_sv_length: "Minimum SV length to consider"
        n_chunks: "Number of chunks to split the VCF into, and subsequently the number of workers"
        n_threads: "Maximum number of threads to use"
        reference_fa: "Reference fasta file"
        skip_solve: "Skip the solve step, only generate input CSV for the solve step"
        tandems_bed: "BED file of tandem repeats"
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
        call hapestry_merge.merge as scattered_merge {
            input:
                vcf_gz = x.left,
                vcf_gz_tbi = x.right,
                bam_not_hardclipped = bam_not_hardclipped,
                interval_max_length = interval_max_length,
                flank_length = flank_length,
                min_sv_length = min_sv_length,
                graphaligner_timeout = graphaligner_timeout,
                min_read_hap_identity = min_read_hap_identity,
                n_threads = n_threads,
                tandems_bed = tandems_bed,
                reference_fa = reference_fa,
                haps_vs_ref_csv = haps_vs_ref_csv,
                force_unique_reads = force_unique_reads,
                skip_solve = skip_solve,
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

    output {
        Array[File] non_sequence_data_tarball = scattered_merge.non_sequence_data_tarball
        Array[File] sequence_data_tarball = scattered_merge.sequence_data_tarball
        File hapestry_vcf = concat_vcfs.concatenated_vcf
        File hapestry_vcf_tbi = concat_vcfs.concatenated_vcf_tbi
    }
}
