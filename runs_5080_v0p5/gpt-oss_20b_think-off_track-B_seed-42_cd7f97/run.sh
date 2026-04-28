#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF=data/ref/chrM.fa
RESULTS=results

# Create results directory if it doesn't exist
mkdir -p "$RESULTS"

# Index reference if not already indexed
if ! [[ -f "${REF}.bwt" && -f "${REF}.pac" && -f "${REF}.ann" && -f "${REF}.amb" && -f "${REF}.sa" ]]; then
    bwa index -p "${REF}" "$REF"
fi
if ! [[ -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

# List of samples
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

# Function to process a single sample
process_sample() {
    local sample="$1"
    local bam="${RESULTS}/${sample}.bam"
    local bai="${bam}.bai"
    local vcf="${RESULTS}/${sample}.vcf.gz"
    local tbi="${vcf}.tbi"

    # Skip if all files exist
    if [[ -f "$bam" && -f "$bai" && -f "$vcf" && -f "$tbi" ]]; then
        return
    fi

    local fq1="data/raw/${sample}_1.fq.gz"
    local fq2="data/raw/${sample}_2.fq.gz"

    # Align reads
    bwa mem -t "$THREADS" "$REF" "$fq1" "$fq2" | \
    samtools view -bS - | \
    samtools sort -@ "$THREADS" -o "$bam"

    # Index BAM
    samtools index "$bam"

    # Variant calling with lofreq
    lofreq call -f "$REF" -o - "$bam" | \
    bcftools view -O z -o "$vcf"

    # Index VCF
    tabix -p vcf "$vcf"
}

# Process each sample
for s in "${SAMPLES[@]}"; do
    process_sample "$s"
done

# Create collapsed table if not exists
COLLAPSED="${RESULTS}/collapsed.tsv"
if [[ ! -f "$COLLAPSED" ]]; then
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED"
    for s in "${SAMPLES[@]}"; do
        bcftools query -f '%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%AF\n' "${RESULTS}/${s}.vcf.gz" >> "$COLLAPSED"
    done
fi