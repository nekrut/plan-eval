#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
mkdir -p results

# 1. Reference indexing — BWA (idempotent)
if [[ ! -f data/ref/chrM.fa.bwt ]]; then
    bwa index data/ref/chrM.fa
fi

# 2. Reference indexing — samtools faidx (idempotent)
if [[ ! -f data/ref/chrM.fa.fai ]]; then
    samtools faidx data/ref/chrM.fa
fi

# 3-6. Per-sample alignment, sorting, indexing, calling, compression
for sample in "${SAMPLES[@]}"; do
    bam_file="results/${sample}.bam"
    vcf_gz_file="results/${sample}.vcf.gz"

    # Skip if BAM and VCF already exist
    if [[ -f "$bam_file" && -f "${bam_file}.bai" && -f "$vcf_gz_file" && -f "${vcf_gz_file}.tbi" ]]; then
        continue
    fi

    # Alignment + sort
    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" data/ref/chrM.fa "data/raw/${sample}_1.fq.gz" "data/raw/${sample}_2.fq.gz" | samtools sort -@ "$THREADS" -o "$bam_file" -

    # BAM index
    samtools index -@ "$THREADS" "$bam_file"

    # Variant calling — LoFreq
    lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o "results/${sample}.vcf" "$bam_file"

    # VCF compression
    bgzip -f "results/${sample}.vcf"

    # Tabix index
    tabix -p vcf "$vcf_gz_file"
done

# 7. Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
for sample in "${SAMPLES[@]}"; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${sample}.vcf.gz" | awk -v s="$sample" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done