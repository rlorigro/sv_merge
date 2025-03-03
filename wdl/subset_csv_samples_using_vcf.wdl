version 1.0


task SubsetSamplesFromCsv {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File csv
    }

    Int disk_size = 1 + ceil(2 * (size(vcf_gz, "GiB")))

    command <<<
        set -euxo pipefail


        # Extract sample names from the VCF and store them in a list
        vcf_samples=$(bcftools query -l ~{vcf_gz})

        # Filter the CSV file to include only those sample names that are in the VCF
        awk -F, 'NR==FNR{samples[$1]; next} $1 in samples' <(echo "$vcf_samples") ~{csv} > output.csv

    >>>

    output {
        File subsampled_csv = "output.csv"
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


workflow SubsetSamplesFromCsvUsingVcf {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File csv
    }

    call SubsetSamplesFromCsv as subset { input:
        vcf_gz = vcf_gz,
        vcf_gz_tbi = vcf_gz_tbi,
        csv = csv
    }

    output {
        File subsampled_csv = subset.subsampled_csv
    }
}
