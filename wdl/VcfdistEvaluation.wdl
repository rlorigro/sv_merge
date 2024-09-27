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
        f.write(f'{mean}\n')
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
        File eval_vcf
        String region
        File reference_fasta
        File reference_fasta_fai

        File vcfdist_bed_file
        String? vcfdist_extra_args
    }

    scatter (sample in samples) {
        # First clean the vcf
        call clean_vcf.CleanVcfAlleles as clean{
            input:
                vcf_gz = eval_vcf,
                ref_fasta = reference_fasta
        }

        call SubsetSampleFromVcf as SubsetSampleFromVcfEval { input:
            vcf = clean.output_vcf_gz,
            sample = sample,
            region = region,
            reference_fasta_fai = reference_fasta_fai
        }

        call SubsetSampleFromVcf as SubsetSampleFromVcfTruth { input:
            vcf = truth_vcf,
            sample = sample,
            region = region,
            reference_fasta_fai = reference_fasta_fai
        }

        call Vcfdist { input:
            sample = sample,
            eval_vcf = SubsetSampleFromVcfEval.single_sample_vcf,
            truth_vcf = SubsetSampleFromVcfTruth.single_sample_vcf,
            bed_file = vcfdist_bed_file,
            reference_fasta = reference_fasta,
            extra_args = vcfdist_extra_args
        }
    }

    scatter (x in Vcfdist.outputs) {
        Float SV_PREC = x.SV_PREC
        Float SV_RECALL = x.SV_RECALL
        Float SV_F1_SCORE = x.SV_F1_SCORE
    }

    call ComputeAverage as avg_p { input:
        x = SV_PREC
    }

    call ComputeAverage as avg_r { input:
        x = SV_RECALL
    }

    call ComputeAverage as avg_f { input:
        x = SV_F1_SCORE
    }

    output {
        # per-sample
        Array[VcfdistOutputs] vcfdist_outputs_per_sample = Vcfdist.outputs

        Float SV_PREC_avg = avg_p.y
        Float SV_RECALL_avg = avg_r.y
        Float SV_F1_SCORE_avg = avg_f.y
    }
}

task SubsetSampleFromVcf {
    input {
        File vcf
        String sample
        String region
        File reference_fasta_fai
    }

    Int disk_size = 1 + ceil(2 * (size(vcf, "GiB")))

    command <<<
        set -euxo pipefail

        bcftools index ~{vcf}
        bcftools view ~{vcf} \
            -s ~{sample} \
            -r ~{region} \
            -Oz -o ~{sample}.subset.g.vcf.gz
        bcftools reheader ~{sample}.subset.g.vcf.gz \
            --fai ~{reference_fasta_fai} \
            -o ~{sample}.subset.reheadered.g.vcf.gz
        bcftools index -t ~{sample}.subset.reheadered.g.vcf.gz
    >>>
    
    output {
        File single_sample_vcf = "~{sample}.subset.reheadered.g.vcf.gz"
        File single_sample_vcf_tbi = "~{sample}.subset.reheadered.g.vcf.gz.tbi"
    }

    runtime {
        cpu: 1
        memory: "64 GiB"
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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

        # Check if input CSV file is provided
        if [ -z "$input_csv" ]; then
            echo "Usage: $0 <input_csv>"
            exit 1
        fi

        # Extract the row where VAR_TYPE = 'SV' and THRESHOLD = 'BEST'
        row=$(awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) header[i]=$i} $1=="SV" && $2=="BEST" {for (i=1; i<=NF; i++) print header[i], $i}' "$input_csv")

        # Extract PREC, RECALL, and F1_SCORE values
        SV_PREC=$(echo "$row" | awk '/PREC/ {print $2}')
        SV_RECALL=$(echo "$row" | awk '/RECALL/ {print $2}')
        SV_F1_SCORE=$(echo "$row" | awk '/F1_SCORE/ {print $2}')

        # Print the extracted values
        echo "SV_PREC: $SV_PREC"
        echo "SV_RECALL: $SV_RECALL"
        echo "SV_F1_SCORE: $SV_F1_SCORE"

        # ^thanks AI

        # write the values to 3 files
        echo "$SV_PREC" > SV_PREC
        echo "$SV_RECALL" > SV_RECALL
        echo "$SV_F1_SCORE" > SV_F1_SCORE

        for tsv in $(ls *.tsv); do mv $tsv ~{sample}.$tsv; done
        mv summary.vcf ~{sample}.summary.vcf
    >>>

    output {
        VcfdistOutputs outputs = {
            "summary_vcf": "~{sample}.summary.vcf",
            "precision_recall_summary_tsv": "~{sample}.precision-recall-summary.tsv",
            "precision_recall_tsv": "~{sample}.precision-recall.tsv",
            "query_tsv": "~{sample}.query.tsv",
            "truth_tsv": "~{sample}.truth.tsv",
            "phasing_summary_tsv": "~{sample}.phasing-summary.tsv",
            "switchflips_tsv": "~{sample}.switchflips.tsv",
            "superclusters_tsv": "~{sample}.superclusters.tsv",
            "phase_blocks_tsv": "~{sample}.phase-blocks.tsv",
            "SV_PREC":read_float("SV_PREC"),
            "SV_RECALL":read_float("SV_RECALL"),
            "SV_F1_SCORE":read_float("SV_F1_SCORE"),
        }
    }

    runtime {
        docker: "timd1/vcfdist:v2.5.3"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}
