version 1.0


task SubsetNSamples {
    input {
        File vcf_gz
        File vcf_gz_tbi
        Int n
    }

    Int disk_size = 1 + ceil(2 * (size(vcf_gz, "GiB")))

    command <<<
        set -euxo pipefail

        bcftools view -Oz -s $(bcftools query -l ~{vcf_gz} | head -n ~{n} | paste -sd, -) ~{vcf_gz} -o first_n_samples.vcf.gz
        bcftools index -t --threads 4 first_n_samples.vcf.gz
    >>>

    output {
        File first_n_samples_vcf_gz = "first_n_samples.vcf.gz"
        File first_n_samples_vcf_gz_tbi = "first_n_samples.vcf.gz.tbi"
    }

    runtime {
        cpu: 4
        memory: "64 GiB"
        disks: "local-disk " + disk_size + " SSD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "staphb/bcftools:1.20"
    }
}


workflow SubsetNSamplesFromVcf {
    input {
        Int n
        File vcf_gz
        File vcf_gz_tbi
    }

    call SubsetNSamples as subset { input:
        vcf_gz = vcf_gz,
        vcf_gz_tbi = vcf_gz_tbi,
        n = n
    }

    output {
        File first_n_samples_vcf_gz = subset.first_n_samples_vcf_gz
        File first_n_samples_vcf_gz_tbi = subset.first_n_samples_vcf_gz_tbi
    }
}
