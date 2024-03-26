#!/usr/bin/env bash

VCF="/home/ryan/data/test_hapestry/vcf/hg002_nist/truvari_8x_bench/truvari_output_confident_only/truvari_HG002_collapsed_NIST_benchnotated_confident_only.vcf.gz"
LABEL="truvari_hg002_benchnotate_max_multiregion"

CSV_READS="/home/ryan/data/test_hapestry/bam/hg002_nist/HG002_8x_hprc_hifi_vs_grch38.csv"
CSV_REF="/home/ryan/data/test_hapestry/bam/hg002_nist/nist_hg002_vs_hg38.csv"
BED="/home/ryan/data/human/reference/human_GRCh38_no_alt_analysis_set.trf.bed"
REF="/home/ryan/data/human/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

MAX_LENGTH=50000

run_dual_annotation(){
  # replace region string spaces with underscores
  REGION_SUFFIX=$(echo ${REGION_STRING} | tr ' ' '_')
  OUT_DIR_READS="/home/ryan/data/test_hapestry/run/${LABEL}_${REGION_SUFFIX}_reads"
  OUT_DIR_REF="/home/ryan/data/test_hapestry/run/${LABEL}_${REGION_SUFFIX}_ref"

  REF_CHR=${REF//.fa*/"_${REGION_SUFFIX}.fasta"}

  VCF_CHR=${VCF//".vcf.gz"/"_${REGION_SUFFIX}.vcf.gz"}
  VCF_NORM=${VCF_CHR//".vcf.gz"/"_norm.vcf"}

  VCF_NAME=$(basename ${VCF_NORM})
  VCF_ANNOTATED=${VCF_NAME//".vcf"/""}
  VCF_ANNOTATED=$(echo ${VCF_ANNOTATED} | tr '.' '_')
  VCF_ANNOTATED=${OUT_DIR_READS}/${VCF_ANNOTATED}"_annotated.vcf"

  echo ${OUT_DIR_READS}
  echo ${OUT_DIR_REF}
  echo ${VCF}
  #echo ${VCF_CHR}
  echo ${VCF_NORM}
  echo ${VCF_NAME}
  echo ${VCF_ANNOTATED}
  #echo ${REF_CHR}

  set -euxo pipefail

  bcftools view ${VCF} ${REGION_STRING} > ${VCF_CHR}
  bcftools norm --multiallelics - ${VCF_CHR} > ${VCF_NORM}

  samtools faidx ${REF} ${REGION_STRING} > ${REF_CHR}

  sync

  /usr/bin/time -v /home/ryan/code/sv_merge/build/annotate \
  --output_dir ${OUT_DIR_READS} \
  --bam_csv ${CSV_READS} \
  --vcf ${VCF_NORM} \
  --tandems ${BED} \
  --ref ${REF_CHR} \
  --interval_max_length ${MAX_LENGTH} \
  --bam_not_hardclipped \
  --flank_length 200 \
  --n_threads 32 \
  --label HAPESTRY_READS

  sync

  /usr/bin/time -v /home/ryan/code/sv_merge/build/annotate \
  --output_dir ${OUT_DIR_REF} \
  --bam_csv ${CSV_REF} \
  --vcf ${VCF_ANNOTATED} \
  --tandems ${BED} \
  --ref ${REF_CHR} \
  --interval_max_length ${MAX_LENGTH} \
  --flank_length 200 \
  --n_threads 32 \
  --label HAPESTRY_REF
}


REGION_STRING="chr1 chr2 chr3"
run_dual_annotation


#region_string="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
#run_dual_annotation

