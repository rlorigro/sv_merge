version 1.0
import "fill_vcf_ref_alleles_from_fasta.wdl" as fill_vcf_ref_alleles


task CleanVcf {
    input {
        File vcf_gz
        File? vcf_tbi
        File ref_fasta
    }

    command {
        # if the tbi is not provided, generate it
        if [ ! -f ${vcf_tbi} ]; then
            bcftools index -t ${vcf_gz}
        fi

        # if the fai does not exist, generate it??
        if [ ! -f ${ref_fasta}.fai ]; then
            samtools faidx ${ref_fasta}
        fi

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


workflow CleanVcfAlleles {
    input {
        File vcf_gz
        File? vcf_tbi
        File ref_fasta
    }

    # First fix the symbolic alleles
    call CleanVcf {
        input:
            vcf_gz = vcf_gz,
            vcf_tbi = vcf_tbi,
            ref_fasta = ref_fasta
    }

    # Finally fix some Ns from the now non-symbolic alleles
    call fill_vcf_ref_alleles.FillFromFasta as fill_from_fasta {
        input:
            vcf_gz = CleanVcf.output_vcf_gz,
            vcf_tbi = CleanVcf.output_vcf_tbi,
            ref_fasta = ref_fasta
    }

    output {
        File output_vcf_gz = fill_from_fasta.output_vcf_gz
        File output_vcf_tbi = fill_from_fasta.output_vcf_tbi
    }
}
