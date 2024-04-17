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
task bcftools_norm {
    input {
        File vcf_gz
        File vcf_gz_tbi

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command <<<
        set -eou pipefail

        if [ -s ~{monitoring_script} ]; then
          bash ~{monitoring_script} > monitoring.log &
        fi

        bcftools norm --multiallelics - -Oz -o temp.vcf.gz ~{vcf_gz}
        mv temp.vcf.gz ~{vcf_gz}

        # Index the vcf
        bcftools index -t ~{vcf_gz}
    >>>

    parameter_meta {
        vcf_gz: "VCF file to annotate"
        vcf_gz_tbi: "Tabix index for the VCF"
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
        File out_vcf_gz = vcf_gz
        File out_vcf_gz_tbi = vcf_gz_tbi
        File? monitoring_log = "monitoring.log"
    }
}


# Define the task
task bcftools_collapse {
    input {
        Array[File] vcfs_gz
        Array[File] vcfs_gz_tbi

        File? bed

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command <<<
        set -eou pipefail

        if [ -s ~{monitoring_script} ]; then
          bash ~{monitoring_script} > monitoring.log &
        fi

        # Use bcftools to concat and collapse the vcfs
        bcftools concat -Oz -o out.vcf.gz --collapse none ~{sep=" " vcfs_gz}

        # Index the vcf
        bcftools index -t out.vcf.gz

        # Use bcftools to subset the vcf by the BED if it was provided
        if ~{defined(bed)}; then
            bcftools view -R ~{bed} -Oz -o temp.vcf.gz out.vcf.gz
            mv temp.vcf.gz out.vcf.gz
            bcftools index -t out.vcf.gz
        fi

    >>>

    parameter_meta {
        bed: "BED file to subset the VCF"
        norm_biallelic: "Normalize the VCF to biallelic sites"
        vcf_gz: "VCF file to annotate"
        vcf_gz_tbi: "Tabix index for the VCF"
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
        File out_vcf_gz = "out.vcf.gz"
        File out_vcf_gz_tbi = "out.vcf.gz.tbi"
        File? monitoring_log = "monitoring.log"
    }
}


workflow bcftools_workflow {
    input {
        Array[File] vcfs_gz
        Array[File] vcfs_gz_tbi
        File? bed
        Boolean norm_biallelic = false

        String docker

        RuntimeAttributes? annotate_runtime_attributes
    }

    parameter_meta {
        bed: "BED file to subset the VCF"
        norm_biallelic: "Normalize the VCF to biallelic sites"
        vcfs_gz: "Array of VCF files to annotate"
        vcfs_gz_tbi: "Tabix indexes for the VCFs"
    }

    call bcftools_collapse as bcf_collapse {
        input:
            vcfs_gz = vcfs_gz,
            vcfs_gz_tbi = vcfs_gz_tbi,
            bed = bed,
            docker = docker
    }

    if (norm_biallelic) {
        call bcftools_norm as bcf_norm {
            input:
                vcf_gz = bcf_collapse.out_vcf_gz,
                vcf_gz_tbi = bcf_collapse.out_vcf_gz_tbi,
                docker = docker
        }
    }

    output {
        File bcftools_vcf_gz = select_first([bcf_norm.out_vcf_gz, bcf_collapse.out_vcf_gz])
        File bcftools_vcf_gz_tbi = select_first([bcf_norm.out_vcf_gz_tbi, bcf_collapse.out_vcf_gz_tbi])
    }
}
