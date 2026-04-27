#!/usr/bin/env bash
set -euo pipefail

THREADS=4
RESULTS_DIR="results"
DATA_RAW="data/raw"
DATA_REF="data/ref"

mkdir -p "$RESULTS_DIR"

# Check if reference is indexed, if not, index it
if [ ! -f "$DATA_REF/chrM.fa.fai" ]; then
    samtools faidx "$DATA_REF/chrM.fa"
fi

# Create a list of samples
samples=()
for file in "$DATA_RAW"/*.fq.gz; do
    sample="${file##*/}"
    sample="${sample%_1.fq.gz}"
    if [[ ! " ${samples[@]} " =~ " $sample " ]]; then
        samples+=("$sample")
    fi
done

# Process each sample
for sample in "${samples[@]}"; do
    # Define output paths
    bam="$RESULTS_DIR/${sample}.bam"
    bai="$RESULTS_DIR/${sample}.bam.bai"
    vcf="$RESULTS_DIR/${sample}.vcf.gz"
    tbi="$RESULTS_DIR/${sample}.vcf.gz.tbi"

    # Skip if all outputs already exist
    if [ -f "$bam" ] && [ -f "$bai" ] && [ -f "$vcf" ] && [ -f "$tbi" ]; then
        continue
    fi

    # Align with bwa
    if [ ! -f "$bam" ]; then
        bwa mem -t "$THREADS" "$DATA_REF/chrM.fa" "$DATA_RAW/${sample}_1.fq.gz" "$DATA_RAW/${sample}_2.fq.gz" | \
        samtools view -Sb - | \
        samtools sort -@ "$THREADS" -o "$bam" -
    fi

    # Index BAM
    if [ ! -f "$bai" ]; then
        samtools index "$bam" "$bai"
    fi

    # Variant calling with lofreq
    if [ ! -f "$vcf" ]; then
        lofreq mpileup -f "$DATA_REF/chrM.fa" "$bam" | \
        lofreq call -f "$DATA_REF/chrM.fa" -o "$vcf" -
    fi

    # Index VCF
    if [ ! -f "$tbi" ]; then
        tabix -p vcf "$vcf"
    fi
done

# Generate collapsed.tsv
collapsed="$RESULTS_DIR/collapsed.tsv"
if [ ! -f "$collapsed" ]; then
    # Extract variants from all VCFs
    for sample in "${samples[@]}"; do
        vcf="$RESULTS_DIR/${sample}.vcf.gz"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' "$vcf" | \
        awk -v sample="$sample" '{print sample, $1, $2, $3, $4, $5}' >> "$collapsed"
    done

    # Sort and deduplicate
    sort -k2,2 -k3,3n "$collapsed" | \
    awk 'seen[$1,$2,$3,$4,$5]++ == 0' > "$collapsed.tmp" && mv "$collapsed.tmp" "$collapsed"
fi