#!/usr/bin/env bash

#VCF="/home/ryan/data/test_hapestry/vcf/hg002_nist/HG002.truvari_collapsed.vcf.gz"
VCF="/home/ryan/data/test_hapestry/vcf/hg002_nist/truvari_8x_bench/truvari_output/truvari_HG002_collapsed_NIST_benchnotated.vcf.gz"
LABEL="truvari_hg002_nist_tandem_max50kbp_benchnotate_50flank_1_skip"


#LABEL="truvari_hg002_nist"
CSV_READS="/home/ryan/data/test_hapestry/bam/hg002_nist/HG002_8x_hprc_hifi_vs_grch38.csv"
CSV_REF="/home/ryan/data/test_hapestry/bam/hg002_nist/nist_hg002_vs_hg38.csv"
BED="/home/ryan/data/human/reference/human_GRCh38_no_alt_analysis_set.trf.bed"
REF="/home/ryan/data/human/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

MAX_LENGTH=50000

run_dual_annotation(){
  OUT_DIR_READS="/home/ryan/data/test_hapestry/run/${LABEL}_${CHR}_reads"
  OUT_DIR_REF="/home/ryan/data/test_hapestry/run/${LABEL}_${CHR}_ref"

  REF_CHR=${REF//.fa*/"_${CHR}.fasta"}

  VCF_CHR=${VCF//".vcf.gz"/"_${CHR}.vcf.gz"}
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

  bcftools view ${VCF} ${CHR} > ${VCF_CHR}
  bcftools norm --multiallelics - ${VCF_CHR} > ${VCF_NORM}

  samtools faidx ${REF} ${CHR} > ${REF_CHR}

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


#for CHR in "chr10"; do
#  run_dual_annotation
#done

for CHR in "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10"; do
  run_dual_annotation
done

#for CHR in "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20"; do
#  run_dual_annotation
#done
