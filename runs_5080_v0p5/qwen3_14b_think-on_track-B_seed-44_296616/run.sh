#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

# Index reference if not already indexed
if [ ! -f data/ref/chrM.fa.fai ]; then
    samtools faidx data/ref/chrM.fa
fi

samples=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

for sample in "${samples[@]}"; do
    # Check if BAM and BAI exist
    if [ ! -f results/${sample}.bam ] || [ ! -f results/${sample}.bam.bai ]; then
        # Align with bwa
        bwa mem data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
        samtools view -Sb - | \
        samtools sort -o results/${sample}.bam -T results/${sample}_tmp
        samtools index results/${sample}.bam
    fi

    # Check if VCF.gz and tbi exist
    if [ ! -f results/${sample}.vcf.gz ] || [ ! -f results/${sample}.vcf.gz.tbi ]; then
        # Call variants with lofreq
        lofreq call -f data/ref/chrM.fa results/${sample}.bam -o results/${sample}.vcf
        bgzip results/${sample}.vcf
        tabix -p vcf results/${sample}.vcf.gz
    fi
done

# Generate collapsed.tsv
if [ ! -f results/collapsed.tsv ]; then
    > results/collapsed.tsv
    for sample in "${samples[@]}"; do
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz | \
        awk -v sample="$sample" '{print sample "\t" $0}' >> results/collapsed.tsv
    done
fi