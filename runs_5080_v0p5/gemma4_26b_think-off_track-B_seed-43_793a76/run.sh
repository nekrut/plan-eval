#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

REF="data/ref/chrM.fa"
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
THREADS=4

if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "$REF"
fi

if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

for SAMPLE in "${SAMPLES[@]}"; do
    BAM="results/${SAMPLE}.bam"
    VCF="results/${SAMPLE}.vcf.gz"

    if [[ -f "$BAM" && -f "$BAM.bai" ]]; then
        continue
    fi

    bwa mem -t "$THREADS" "$REF" \
        "data/raw/${SAMPLE}_1.fq.gz" \
        "data/raw/${SAMPLE}_2.fqgs" 2>/dev/null || \
    bwa mem -t "$THREADS" "$REF" \
        "data/raw/${SAMPLE}_1.fq.gz" \
        "data/raw/${SAMPLE}_2.fq.gz" | \
    samtools view -Usa0 - | \
    samtools sort -@ "$THREADS" -o "$BAM"

    samtools index "$BAM"
done

for SAMPLE in "${SAMPLES[@]}"; do
    VCF="results/${SAMPLE}.vcf.gz"
    BAM="results/${SAMPLE}.bam"

    if [[ -f "$VCF" && -f "$VCF.tbi" ]]; then
        continue
    fi

    lofreq call -f "$REF" -o "$VCF" "$BAM"
    
    # lofreq output might not be compressed or indexed correctly for downstream
    if [[ ! -f "$VCF.tbi" ]]; then
        bgzip -c "$VCF" > "${VCF}.tmp"
        mv "${VCF}.tmp" "$VCF"
        tabix -p vcf "$VCF"
    fi
done

# Create collapsed table
COLLAPSED="results/collapsed.tsv"
if [[ ! -f "$COLLAPSED" ]]; then
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED"
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF="results/${SAMPLE}.vcf.gz"
        # Use bcftools query to extract info. 
        # lofreq VCFs contain AF in INFO field.
        bcftools query -f "$SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF" >> "$COLLAPSED"
    done
fi

exit 0