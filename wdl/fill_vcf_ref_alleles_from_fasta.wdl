version 1.0

workflow FillVcfRefAlleles {
    input {
        File vcf_gz
        File vcf_tbi
        File ref_fasta
    }

    call FillFromFasta {
        input:
            vcf_gz = vcf_gz,
            vcf_tbi = vcf_tbi,
            ref_fasta = ref_fasta
    }

    output {
        File output_vcf_gz = FillFromFasta.output_vcf_gz
        File output_vcf_tbi = FillFromFasta.output_vcf_tbi
    }
}

task FillFromFasta {
    input {
        File vcf_gz
        File vcf_tbi
        File ref_fasta
    }

    command {
        bcftools +fill-from-fasta -Oz ${vcf_gz} -- -c REF -f ${ref_fasta} > output.vcf.gz
        bcftools index -t output.vcf.gz
    }

    output {
        File output_vcf_gz = "output.vcf.gz"
        File output_vcf_tbi = "output.vcf.gz.tbi"
    }

    runtime {
        docker: "staphpb/bcftools:latest"
        memory: "4G"
        cpu: 1
    }
}