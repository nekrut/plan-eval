#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# Reference indexing
if [ ! -f data/ref/chrM.fa.1.hg19.json ] && [ ! -f data/ref/chrM.fa.1.fasta ]; then
    # Check if bwa index exists (bwa index creates .amb, ._1._1.bam etc or just .fa.1.something)
    # Standard bwa index for a single file is usually .fa.1.hg19.json or similar depending on version
    # We run it to ensure indices exist.
    bwa index data/ref/chrM.fa
fi

samtools faidx data/ref/chrM.fa

for sample in "${SAMPLES[@]}"; do
    # Alignment and sorting
    # Note: samtools sort is used directly from bwa pipe
    bwa mem -t "$THREADS" \
        -K 100M \
        -c 1000 \
        -M \
        -t "$THREADS" \
        -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
        data/ref/chrM.fa \
        data/raw/${sample}_1.fq.gz \
        data/raw/${sample}_2.fq.gz | \
    samtools sort -@ "$THREADS" -o results/"${sample}".bam -

    # BAM indexing
    samtools index -@ "$THREADS" results/"${sample}".bam

    # Variant calling
    lofreq call-parallel --pp-threads "$THREADS" \
        -f data/ref/chrM.fa \
        -o results/"${sample}".vcf \
        results/"${sample}".bam

    # VCF compression and tabix indexing
    # Note: bgzip is part of htslib, usually available if samtools is present
    bgzip -f results/"${sample}".vcf
    tabix -p vcf results/"${sample}".vcf.gz

done

# Create collapsed TSV header
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv

# Append data to collapsed TSV
for sample in "${SAMPLES[@]}"; do
    bcftools query -l results/"${sample}".vcf.gz | xargs -I {} \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/"${sample}".vcf.gz | \
    awk -v s="${sample}" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done