version 1.0
import "clean_vcf.wdl" as clean_vcf


struct VcfdistOutputs {
    File summary_vcf
    File precision_recall_summary_tsv
    File precision_recall_tsv
    File query_tsv
    File truth_tsv
    File phasing_summary_tsv
    File switchflips_tsv
    File superclusters_tsv
    File phase_blocks_tsv
    Float SV_PREC
    Float SV_RECALL
    Float SV_F1_SCORE
}


# task to average any Array[Float]
task ComputeAverage {
    input {
        Array[Float] x
    }

    command <<<
python <<CODE
arr = ["~{sep="\", \"" x}"]
floats = [float(x) for x in arr]
mean = sum(floats) / len(floats)
with open('mean.txt', 'w') as file:
    file.write(f'{mean}\n')
CODE
    >>>

    runtime {
        cpu: 1
        memory: "8 GiB"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "python:3.9.20-slim-bullseye"
    }

    output {
        Float y = read_float("mean.txt")
    }
}


workflow VcfdistEvaluation {
    input {
        Array[String] samples
        File truth_vcf
        File? truth_vcf_tbi
        File eval_vcf
        File? eval_vcf_tbi
        String region
        File reference_fasta
        File reference_fasta_fai

        File vcfdist_bed_file
        String? vcfdist_extra_args
    }

    # First clean the vcf
    call clean_vcf.CleanVcfAlleles as clean{
        input:
            vcf_gz = eval_vcf,
            ref_fasta = reference_fasta
    }

    scatter (sample in samples) {

        call SubsetSampleFromVcf as SubsetSampleFromVcfEval { input:
            vcf = clean.output_vcf_gz,
            vcf_tbi = clean.output_vcf_tbi,
            sample = sample,
            region = region,
            reference_fasta_fai = reference_fasta_fai
        }

        call SubsetSampleFromVcf as SubsetSampleFromVcfTruth { input:
            vcf = truth_vcf,
            vcf_tbi = truth_vcf_tbi,
            sample = sample,
            region = region,
            reference_fasta_fai = reference_fasta_fai
        }

        call Vcfdist as vcfdist { input:
            sample = sample,
            eval_vcf = SubsetSampleFromVcfEval.single_sample_vcf,
            truth_vcf = SubsetSampleFromVcfTruth.single_sample_vcf,
            bed_file = vcfdist_bed_file,
            reference_fasta = reference_fasta,
            extra_args = vcfdist_extra_args
        }
    }

    call ComputeAverage as avg_p { input:
        x = vcfdist.SV_PREC
    }

    call ComputeAverage as avg_r { input:
        x = vcfdist.SV_RECALL
    }

    call ComputeAverage as avg_f { input:
        x = vcfdist.SV_F1_SCORE
    }

    output {
        # per-sample
        Array[File] vcfdist_summary_vcf = vcfdist.summary_vcf
        Array[File] vcfdist_precision_recall_summary_tsv = vcfdist.precision_recall_summary_tsv
        Array[File] vcfdist_superclusters_tsv = vcfdist.superclusters_tsv

        Float SV_PREC_avg = avg_p.y
        Float SV_RECALL_avg = avg_r.y
        Float SV_F1_SCORE_avg = avg_f.y
    }
}

task SubsetSampleFromVcf {
    input {
        File vcf
        File? vcf_tbi
        String sample
        String region
        File reference_fasta_fai
    }

    Int disk_size = 1 + ceil(2 * (size(vcf, "GiB")))

    command <<<
        set -euxo pipefail

        if ~{!defined(vcf_tbi)}; then
            bcftools index --threads 4 ~{vcf}
        fi

        bcftools view ~{vcf} \
            -r ~{region} \
            --threads 4 \
            -Oz -o region.vcf.gz

        bcftools index -t --threads 4 region.vcf.gz

        bcftools view region.vcf.gz \
            --threads 4 \
            -s ~{sample} \
            -r ~{region} \
            -Oz -o ~{sample}.subset.g.vcf.gz

        bcftools reheader ~{sample}.subset.g.vcf.gz \
            --fai ~{reference_fasta_fai} \
            -o ~{sample}.subset.reheadered.g.vcf.gz

        bcftools index --threads 4 -t ~{sample}.subset.reheadered.g.vcf.gz
    >>>

    output {
        File single_sample_vcf = "~{sample}.subset.reheadered.g.vcf.gz"
        File single_sample_vcf_tbi = "~{sample}.subset.reheadered.g.vcf.gz.tbi"
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

task Vcfdist {
    input {
        String sample
        File eval_vcf
        File truth_vcf
        File bed_file
        File reference_fasta
        String? extra_args
        Int verbosity = 1

        Int disk_size_gb = ceil(size(truth_vcf, "GiB") + 10)
        Int mem_gb = 16
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        vcfdist \
            ~{eval_vcf} \
            ~{truth_vcf} \
            ~{reference_fasta} \
            -b ~{bed_file} \
            -v ~{verbosity} \
            ~{extra_args}

        # Extract the row where VAR_TYPE = 'SV' and THRESHOLD = 'BEST'
        row=$(awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) header[i]=$i} $1=="SV" && $2=="BEST" {for (i=1; i<=NF; i++) print header[i], $i}' precision-recall-summary.tsv)

        # Extract PREC, RECALL, and F1_SCORE values
        SV_PREC=$(echo "$row" | awk '/PREC/ {print $2}')
        SV_RECALL=$(echo "$row" | awk '/RECALL/ {print $2}')
        SV_F1_SCORE=$(echo "$row" | awk '/F1_SCORE/ {print $2}')

        # Print the extracted values
        echo "SV_PREC: $SV_PREC"
        echo "SV_RECALL: $SV_RECALL"
        echo "SV_F1_SCORE: $SV_F1_SCORE"

        # write the values to 3 files
        echo "$SV_PREC" > SV_PREC
        echo "$SV_RECALL" > SV_RECALL
        echo "$SV_F1_SCORE" > SV_F1_SCORE

        for tsv in $(ls *.tsv); do mv $tsv ~{sample}.$tsv; done
        mv summary.vcf ~{sample}.summary.vcf
    >>>

    output {
        File summary_vcf = "~{sample}.summary.vcf"
        File precision_recall_summary_tsv = "~{sample}.precision-recall-summary.tsv"
        File precision_recall_tsv = "~{sample}.precision-recall.tsv"
        File query_tsv = "~{sample}.query.tsv"
        File truth_tsv = "~{sample}.truth.tsv"
        File phasing_summary_tsv = "~{sample}.phasing-summary.tsv"
        File switchflips_tsv = "~{sample}.switchflips.tsv"
        File superclusters_tsv = "~{sample}.superclusters.tsv"
        File phase_blocks_tsv = "~{sample}.phase-blocks.tsv"
        Float SV_PREC =read_float("SV_PREC")
        Float SV_RECALL =read_float("SV_RECALL")
        Float SV_F1_SCORE =read_float("SV_F1_SCORE")
    }

    runtime {
        docker: "timd1/vcfdist:v2.5.3"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}
