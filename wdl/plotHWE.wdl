version 1.0


task subset_vcf {
    input {
        File vcf_gz
        File vcf_tbi
        String? region
    }

    command {
        if [[ ~{defined(region)} == true ]]; then
        echo "Subsetting to region ~{region}"
        bcftools view -r ~{region} ~{vcf_gz} -Oz -o subset.vcf.gz
        tabix -p vcf subset.vcf.gz
        else
        echo "No region specified. Copying input VCF."
        cp ~{vcf_gz} subset.vcf.gz
        cp ~{vcf_tbi} subset.vcf.gz.tbi
        fi
    }

    output {
        File vcf_out = "subset.vcf.gz"
        File tbi_out = "subset.vcf.gz.tbi"
    }

    runtime {
        docker: "staphb/bcftools:1.21"
        memory: "4G"
        cpu: 1
    }
}


task get_gt_counts {
    input {
        File vcf_gz
    }

    command {
        set -euo pipefail

        # Download and compile the Java script
        curl -L -o PlotHW.java https://raw.githubusercontent.com/rlorigro/sv_merge/dev/wdl/PlotHW.java
        javac PlotHW.java

        mkdir -p output
        java -cp . PlotHW ${vcf_gz} output
    }

    output {
        File output_txt = "output/genotypes_all.txt"
    }

    runtime {
        docker: "openjdk:25-jdk-bullseye"
        memory: "4G"
        cpu: 1
    }
}


task plot_hwe_from_counts {
    input {
        File input_matrix
    }

    command {
        set -euo pipefail

        # Download the R script
        wget https://raw.githubusercontent.com/rlorigro/sv_merge/dev/wdl/PlotHW.r

        Rscript PlotHW.r ${input_matrix} hwe_plot.png
    }

    output {
        File output_png = "hwe_plot.png"
    }

    runtime {
        docker: "r-base:4.4.2"
        memory: "4G"
        cpu: 1
    }
}


workflow plot_hwe {
    input {
        File input_vcf_gz
        File input_vcf_tbi
        String? region
    }

    call subset_vcf {
        input:
            vcf_gz = input_vcf_gz,
            vcf_tbi = input_vcf_tbi,
            region = region
    }

    call get_gt_counts {
        input:
            vcf_gz = subset_vcf.vcf_out
    }

    call plot_hwe_from_counts {
        input:
            input_matrix = get_gt_counts.output_txt
    }

    output {
        File hwe_plot = plot_hwe_from_counts.output_png
    }
}
