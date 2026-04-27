#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# Reference indexing
[[ -f data/ref/chrM.fa.bwt ]] || bwa index data/ref/chrM.fa
[[ -f data/ref/chrM.fa.fai ]] || samtools faidx data/ref/chrM.fa

for sample in "${SAMPLES[@]}"; do
    # Alignment and sort
    [[ -f results/${sample}.bam ]] || {
        bwa mem -t ${THREADS} -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz \
            | samtools sort -@ ${THREADS} -o results/${sample}.bam -
    }

    # BAM index
    [[ -f results/${sample}.bam.bai ]] || samtools index -@ ${THREADS} results/${sample}.bam

    # Variant calling
    [[ -f results/${sample}.vcf.gz.tbi ]] || {
        lofreq call-parallel --pp-threads ${THREADS} -f data/ref/chrM.fa -o results/${sample}.vcf results/${sample}.bam
    }

    # Compression and tabix index
    [[ -f results/${sample}.vcf.gz.tbi ]] || {
        samtools bgzip -f results/${sample}.vcf
        tabix -p vcf results/${sample}.vcf.gz
    }
done

# Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
for sample in "${SAMPLES[@]}"; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz | \
    while IFS=$'\t' read -r chrom pos ref alt af; do
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$chrom" "$pos" "$ref" "$alt" "$af"
    done >> results/collapsed.tsv
done