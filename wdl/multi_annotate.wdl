version 1.0


import "https://raw.githubusercontent.com/fabio-cunial/callset_integration/7027bf4f282ab902f306c3fe2d086db18fcd7883/wdl/GetRegenotypedVcfSniffles.wdl" as sniffles
import "https://github.com/fabio-cunial/callset_integration/raw/4c3a061b5db1eb99b08dae6564d75a76ad87e933/wdl/AddTruvariAnnotations.wdl" as truvari
import "hapestry_annotate.wdl" as hapestry_annotate
import "https://raw.githubusercontent.com/fabio-cunial/callset_integration/main/wdl/RepeatAnnotation.wdl" as repeat_annotation


workflow multiannotate {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File confident_bed

        String sample_id
        Int max_length = 50000
        Int min_length = 0
        Int flank_length = 200
        Int n_threads

        # annotations
        File tandems_bed
        File segdup_bed
        File telomere_bed
        File centromere_bed

        File reference_fa
        File reference_fai
        File truth_vcf_gz
        File truth_vcf_gz_tbi
        File haps_vs_ref_csv
        File reads_vs_ref_csv
        File alignments_bam
        File alignments_bai
        Boolean force_unique_reads = false
        String truvari_bench_flags = "--sizemin 0 --sizefilt 0 --sizemax 1000000 --pctsize 0.9 --pctseq 0.9 --pick multi"

        String docker

        RuntimeAttributes? annotate_runtime_attributes
    }

    parameter_meta {
        min_score: "Minimum identity observed spanning a variant to consider a true positive w.r.t. assembly alignments"
        interval_max_length: "Maximum length of each window evaluated"
        flank_length: "Length of flanking sequence to include in each window"
        n_threads: "Maximum number of threads to use"
        tandems_bed: "BED file of tandem repeats"
        segdup_bed: "BED file of segmental duplications"
        telomere_bed: "BED file of telomeres"
        centromere_bed: "BED file of centromeres"
        reference_fa: "Reference fasta file"
        haps_vs_ref_csv: "CSV file of haplotype vs reference BAMs"
        force_unique_reads: "Force unique aligned sequence names among multiple BAMs to prevent collisions"
        bam_not_hardclipped: "If the bam is GUARANTEED not to contain any hardclips, use this flag to trigger much simpler/faster fetching process"
        annotation_label: "Name to give the INFO field in the VCF for annotations, usually upper case"
    }

    call hapestry_annotate.annotate as annotate_reads {
        input:
            vcf_gz = vcf_gz,
            vcf_gz_tbi = vcf_gz_tbi,
            confident_bed = confident_bed,
            bam_not_hardclipped = true,
            interval_max_length = max_length,
            flank_length = flank_length,
            n_threads = n_threads,
            tandems_bed = tandems_bed,
            reference_fa = reference_fa,
            haps_vs_ref_csv = reads_vs_ref_csv,
            force_unique_reads = force_unique_reads,
            annotation_label = "HAPESTRY_READS",
            docker = docker,
            runtime_attributes = annotate_runtime_attributes
    }

    call hapestry_annotate.annotate as annotate_ref {
        input:
            vcf_gz = annotate_reads.annotated_vcf_gz,
            vcf_gz_tbi = annotate_reads.annotated_vcf_gz_tbi,
            confident_bed = confident_bed,
            bam_not_hardclipped = false,
            interval_max_length = max_length,
            flank_length = flank_length,
            n_threads = n_threads,
            tandems_bed = tandems_bed,
            reference_fa = reference_fa,
            haps_vs_ref_csv = haps_vs_ref_csv,
            force_unique_reads = force_unique_reads,
            annotation_label = "HAPESTRY_REF",
            docker = docker,
            runtime_attributes = annotate_runtime_attributes
    }

    call truvari.AddTruvariAnnotationsImpl as truvari_benchnotate {
        input:
            sample_id = sample_id,
            input_vcf_gz = annotate_ref.annotated_vcf_gz,
            input_vcf_gz_tbi = annotate_ref.annotated_vcf_gz_tbi,
            truth_vcf_gz = truth_vcf_gz,
            truth_tbi = truth_vcf_gz_tbi,
            truvari_bench_flags = truvari_bench_flags
    }

    call sniffles.GetRegenotypedVcfImpl as sniffles_annotate {
        input:
            truvari_collapsed_vcf_gz = truvari_benchnotate.annotated_vcf_gz,
            truvari_collapsed_tbi = truvari_benchnotate.annotated_tbi,
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            svlen_min = min_length,
            svlen_max = max_length
    }

    call repeat_annotation.RepeatAnnotationImpl as repeat_annotate {
        input:
            vcf_gz = sniffles_annotate.regenotyped_sniffles,
            vcf_gz_tbi = sniffles_annotate.regenotyped_sniffles_tbi,
            segdup_bed = segdup_bed,
            tr_bed = telomere_bed,
            centromere_bed = centromere_bed
    }

    output {
        File annotated_vcf_gz = repeat_annotate.annotated_vcf
        File annotated_vcf_gz_tbi = repeat_annotate.annotated_tbi
    }
}
