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


task predict_variant_scores {
    input {
        File vcf_gz
        File vcf_tbi
        File model
        Float min_score = 0.15
        Int allow_all_above_length = 2048
        String features_name = "hapestry"

        RuntimeAttributes runtime_attributes = {}
    }

    String repository_dir = "/hapestry/sv_merge"

    command <<<
        git --no-pager --git-dir ~{repository_dir}/.git log --decorate=short --pretty=oneline | head -n 1

        set -euxo pipefail

        # get the directory of the vcf
        vcf_dir=$(dirname ~{vcf_gz})

        # get the tbi directory
        tbi_dir=$(dirname ~{vcf_tbi})

        # if the tbi_dir and the vcf_dir are not the same then move the tbi into the vcf directory
        if [ "$vcf_dir" != "$tbi_dir" ]; then
            mv ~{vcf_tbi} ${vcf_dir}
        fi

        python3 ~{repository_dir}/scripts/predict_variant_scores.py \
        --vcf_paths ~{vcf_gz} \
        --features_name ~{features_name} \
        --model_path ~{model} \
        --output_path output/hapestry_filter_scored.vcf

        bcftools view --threads 4 -Oz output/hapestry_filter_scored.vcf -o output/hapestry_filter_scored.vcf.gz
        bcftools index --threads 4 -t output/hapestry_filter_scored.vcf.gz

        bcftools filter -Oz -i "abs(SVLEN) >= ~{allow_all_above_length} | HAPESTRY_SCORE >= ~{min_score}" output/hapestry_filter_scored.vcf.gz -o output/hapestry_filtered.vcf.gz
        bcftools index --threads 4 -t output/hapestry_filtered.vcf.gz
    >>>

    output {
        File hapestry_filtered_vcf_gz = "output/hapestry_filtered.vcf.gz"
        File hapestry_filtered_vcf_tbi = "output/hapestry_filtered.vcf.gz.tbi"
    }

    runtime {
        docker: "fcunial/hapestry:filter"
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    parameter_meta{
        vcf_gz: "Path to the VCF file"
        vcf_tbi: "Path to the VCF index"
        features_name: "which feature set to use for inference (e.g. 'hapestry' or 'sniffles')"
        model: "Path to the filter model, typically '.pt' file extension"
        min_score: "Minimum score (emitted by Hapestry NN model) of variant to allow in the filtered VCF. Ideally based on a ROC curve."
        allow_all_above_length: "Do NOT filter by min_score for variants above this length. They will still have a HAPESTRY_SCORE INFO field."
    }
}


workflow hapestry_filter {
    input {
        File vcf_gz
        File vcf_tbi
        File model
    }

    call predict_variant_scores {
        input:
            vcf_gz = vcf_gz,
            vcf_tbi = vcf_tbi,
            model = model
    }

    parameter_meta{
        vcf_gz: "Path to the VCF file"
        vcf_tbi: "Path to the VCF index"
        model: "Path to the filter model, typically '.pt' file extension"
    }

    output{
        File hapestry_filtered_vcf_gz = predict_variant_scores.hapestry_filtered_vcf_gz
        File hapestry_filtered_vcf_tbi = predict_variant_scores.hapestry_filtered_vcf_tbi
    }
}
