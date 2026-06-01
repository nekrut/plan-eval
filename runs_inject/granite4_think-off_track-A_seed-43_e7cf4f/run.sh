#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
OUTPUT_DIR="results"
mkdir -p "$OUTPUT_DIR"

# 1. Index reference (bwa & samtools)
if [[ ! -f data/ref/chrM.fa.bwt ]]; then
    bwa index data/ref/chrM.fa
fi

if [[ ! -f data/ref/chrM.fa.fai ]]; then
    samtools faidx data/ref/chrM.fa
fi

# 2. Alignment & sorting per sample (bwa mem + samtools sort)
for sample in "${SAMPLES[@]}"; do
    bwa mem -t "$THREADS" \
        -R "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA" \
        data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz |
        samtools sort -@ "$THREADS" -o "$OUTPUT_DIR/${sample}.bam"
done

# 3. Index BAM
for sample in "${SAMPLES[@]}"; do
    if [[ ! -f "$OUTPUT_DIR/${sample}.bam.bai" ]]; then
        samtools index -@ "$THREADS" "$OUTPUT_DIR/${sample}.bam"
    fi
done

# 4. Variant calling (lofreq)
for sample in "${SAMPLES[@]}"; do
    lofreq call-parallel --pp-threads "$THREADS" \
        -f data/ref/chrM.fa -o "$OUTPUT_DIR/${sample}.vcf" "$OUTPUT_DIR/${sample}.bam"
done

# 5. VCF compression + tabix
for sample in "${SAMPLES[@]}"; do
    if [[ ! -f "$OUTPUT_DIR/${sample}.vcf.gz.tbi" ]]; then
        bgzip -f "$OUTPUT_DIR/${sample}.vcf"
        tabix -p vcf "$OUTPUT_DIR/${sample}.vcf.gz"
    fi
done

# 6. Collapsed TSV (always rebuild)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUTPUT_DIR/collapsed.tsv"

for sample in "${SAMPLES[@]}"; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUTPUT_DIR/${sample}.vcf.gz" |
        awk -v s="$sample" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUTPUT_DIR/collapsed.tsv"
done