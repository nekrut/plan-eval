#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
RESULTS_DIR="results"
RAW_DIR="data/raw"
REF_DIR="data/ref"
REF_FA="$REF_DIR/chrM.fa"

mkdir -p "$RESULTS_DIR"

# 2. Reference indexing (once)
if [[ ! -f "$REF_FA.fai" ]]; then
    samtools faidx "$REF_FA"
fi
if [[ ! -f "$REF_FA.bwt" ]]; then
    bwa index "$REF_FA"
fi

# 3. Per-sample alignment with bwa mem
for sample in "${SAMPLES[@]}"; do
    RG="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
    bwa mem -t "$THREADS" -R "$RG" "$REF_FA" "$RAW_DIR/${sample}_1.fq.gz" "$RAW_DIR/${sample}_2.fq.gz" | \
    samtools sort -@ "$THREADS" -o "$RESULTS_DIR/${sample}.bam"
    samtools index -@ "$THREADS" "$RESULTS_DIR/${sample}.bam"
done

# 6. Variant calling with lofreq call-parallel
for sample in "${SAMPLES[@]}"; do
    lofreq call-parallel --pp-threads "$THREADS" -f "$REF_FA" -o "$RESULTS_DIR/${sample}.vcf" "$RESULTS_DIR/${sample}.bam"
done

# 7. VCF compression and indexing
for sample in "${SAMPLES[@]}"; do
    bgzip -c "$RESULTS_DIR/${sample}.vcf" > "$RESULTS_DIR/${sample}.vcf.gz"
    tabix -p vcf "$RESULTS_DIR/${sample}.vcf.gz"
    rm "$RESULTS_DIR/${sample}.vcf"
done

# 8. Collapse step
HEADER="sample\tchrom\tpos\tref\talt\taf"
{
    echo "$HEADER"
    for sample in "${SAMPLES[@]}"; do
        bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$RESULTS_DIR/${sample}.vcf.gz"
    done
} > "$RESULTS_DIR/collapsed.tsv"