version 1.0

workflow CleanVcfAlleles {
    input {
        File vcf_gz
        File vcf_tbi
        File ref_fasta
    }

    call CleanVcf {
        input:
            vcf_gz = vcf_gz,
            vcf_tbi = vcf_tbi,
            ref_fasta = ref_fasta
    }

    output {
        File output_vcf_gz = CleanVcf.output_vcf_gz
        File output_vcf_tbi = CleanVcf.output_vcf_tbi
    }
}

task CleanVcf {
    input {
        File vcf_gz
        File vcf_tbi
        File ref_fasta
    }

    command {
        # Does the following:
        # - Fixing symbolic ALTs
        # - Uppercasing REF and ALT
        # - Fixing SVLEN
        # - Fixing GTs

        python3 /hgsvc2/resolve_light.py ${vcf_gz} ${ref_fasta} > output.vcf
        bcftools view output.vcf -Oz -o output.vcf.gz
        bcftools index -t output.vcf.gz
    }

    output {
        File output_vcf_gz = "output.vcf.gz"
        File output_vcf_tbi = "output.vcf.gz.tbi"
    }

    runtime {
        docker: "fcunial/callset_integration:latest"
        memory: "4G"
        cpu: 1
    }
}