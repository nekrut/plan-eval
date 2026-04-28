#!/usr/bin/env bash
set -euo pipefail

# Index reference if not already indexed
if [ ! -f data/ref/chrM.fa.bai ]; then
    samtools index data/ref/chrM.fa
fi

samples=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Align and sort BAM
for sample in "${samples[@]}"; do
    if [ -f results/${sample}.bam ]; then
        echo "Skipping alignment for ${sample}"
        continue
    fi
    bwa mem data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -o results/${sample}.bam
    samtools index results/${sample}.bam
done

# Call variants with lofreq
for sample in "${samples[@]}"; do
    if [ -f results/${sample}.vcf ]; then
        echo "Skipping lofreq for ${sample}"
        continue
    fi
    lofreq call -f data/ref/chrM.fa -o results/${sample}.vcf results/${sample}.bam
done

# Convert VCF to gzipped and index
for sample in "${samples[@]}"; do
    if [ -f results/${sample}.vcf.gz ]; then
        echo "Skipping VCF conversion for ${sample}"
        continue
    fi
    bcftools view results/${sample}.vcf -O z -o results/${sample}.vcf.gz
    tabix -p vcf results/${sample}.vcf.gz
done

# Process VCF to get variants and calculate AF
if [ ! -f results/collapsed.tsv ]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv
fi

for sample in "${samples[@]}"; do
    if [ ! -f results/${sample}.vcf.gz ]; then
        echo "VCF file not found for ${sample}, skipping"
        continue
    fi
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AD\t%DP\n' results/${sample}.vcf.gz | awk -v sample=${sample} '{ if (NR > 1) { af = $5 / $6; print sample "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" af } }' >> results/collapsed.tsv
done