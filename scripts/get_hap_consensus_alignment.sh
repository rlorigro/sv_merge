#!/bin/bash

set -euxo

# Define variables for the input VCF, reference fasta, output directory, and sample name
VCF_INPUT=$1
REFERENCE_FASTA=$
OUTPUT_DIR=$3
REGION=$4
SAMPLE_NAME=$5
N_THREADS=$6

# Check if the output directory exists, and create it if not
if [ ! -d "$OUTPUT_DIR" ]; then
  echo "Directory $OUTPUT_DIR does not exist. Creating it now."
  mkdir -p "$OUTPUT_DIR"
else
  echo "Directory $OUTPUT_DIR already exists."
fi


# Process the first VCF file
bcftools view -Oz $VCF_INPUT > $OUTPUT_DIR/merged.vcf.gz
bcftools index -t $OUTPUT_DIR/merged.vcf.gz

# Resolve variants with Python script
python3 /home/ryan/code/sv_merge/scripts/resolve_light.py \
  $OUTPUT_DIR/merged.vcf.gz \
  $REFERENCE_FASTA \
  > $OUTPUT_DIR/merged_resolved.vcf

# Process the resolved VCF
bcftools view -Oz $OUTPUT_DIR/merged_resolved.vcf > $OUTPUT_DIR/merged_resolved.vcf.gz
bcftools index -t $OUTPUT_DIR/merged_resolved.vcf.gz

# Fill VCF from fasta
bcftools +fill-from-fasta -Oz $OUTPUT_DIR/merged_resolved.vcf.gz -- -c REF -f $REFERENCE_FASTA > $OUTPUT_DIR/merged_resolved_filled.vcf.gz
bcftools index -t $OUTPUT_DIR/merged_resolved_filled.vcf.gz

# Generate consensus for haplotype 1
bcftools consensus \
  -s $SAMPLE_NAME \
  -H 1 \
  -f $REFERENCE_FASTA \
  $OUTPUT_DIR/merged_resolved_filled.vcf.gz > $OUTPUT_DIR/${SAMPLE_NAME}_hap_0_full.fasta

# Generate consensus for haplotype 2
bcftools consensus \
  -s $SAMPLE_NAME \
  -H 2 \
  -f $REFERENCE_FASTA \
  $OUTPUT_DIR/merged_resolved_filled.vcf.gz > $OUTPUT_DIR/${SAMPLE_NAME}_hap_1_full.fasta

samtools faidx $OUTPUT_DIR/${SAMPLE_NAME}_hap_0_full.fasta ${REGION} > $OUTPUT_DIR/${SAMPLE_NAME}_hap_0.fasta
samtools faidx $OUTPUT_DIR/${SAMPLE_NAME}_hap_1_full.fasta ${REGION} > $OUTPUT_DIR/${SAMPLE_NAME}_hap_1.fasta

# Run minimap2 to align the haplotype FASTA files to the target reference genome
minimap2 -t ${N_THREADS} -x asm5 -a -L ${REFERENCE_FASTA} $OUTPUT_DIR/${SAMPLE_NAME}_hap_0.fasta $OUTPUT_DIR/${SAMPLE_NAME}_hap_1.fasta > $OUTPUT_DIR/${SAMPLE_NAME}_haps_vs_ref.sam

# Sort the SAM file, convert to BAM, and index it
samtools sort $OUTPUT_DIR/${SAMPLE_NAME}_haps_vs_ref.sam -o $OUTPUT_DIR/${SAMPLE_NAME}_haps_vs_ref.sorted.bam
samtools index $OUTPUT_DIR/${SAMPLE_NAME}_haps_vs_ref.sorted.bam

# Clean up intermediate files
rm $OUTPUT_DIR/merged_resolved.vcf
rm $OUTPUT_DIR/merged_resolved.vcf.gz
rm $OUTPUT_DIR/merged_resolved.vcf.gz.tbi
rm $OUTPUT_DIR/merged_resolved_filled.vcf.gz
rm $OUTPUT_DIR/merged_resolved_filled.vcf.gz.tbi
rm $OUTPUT_DIR/${SAMPLE_NAME}_hap_1_full.fasta
rm $OUTPUT_DIR/${SAMPLE_NAME}_hap_1_full.fasta.fai
rm $OUTPUT_DIR/${SAMPLE_NAME}_hap_0_full.fasta
rm $OUTPUT_DIR/${SAMPLE_NAME}_hap_0_full.fasta.fai
rm $OUTPUT_DIR/${SAMPLE_NAME}_haps_vs_ref.sam

echo "Done"
