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


task MergeBams {
    input {
        Array[File] bams
        Int n_threads
        String output_prefix
        String docker = "staphb/samtools:1.20"

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        samtools merge -@ ~{n_threads} -o ~{output_prefix}.bam ~{sep=" " bams}
        samtools index -@ ~{n_threads} ~{output_prefix}.bam
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
        File merged_bam = output_prefix + ".bam"
        File merged_bam_index = output_prefix + ".bam.bai"
    }
}


workflow MergeWorkflow {
    input {
        Array[File] input_bams
        Int n_threads
        String output_prefix

        RuntimeAttributes? runtime_attributes
    }

    parameter_meta {
        input_bams: "Bams to be merged"
        output_prefix: "Prefix for the merged bam filename"
        n_threads: "Number of threads to use for merging"
    }

    call MergeBams {
        input:
            bams = input_bams,
            n_threads = n_threads,
            output_prefix = output_prefix
    }

    output {
        File merged_bam = MergeBams.merged_bam
        File merged_bam_index = MergeBams.merged_bam_index
    }
}
