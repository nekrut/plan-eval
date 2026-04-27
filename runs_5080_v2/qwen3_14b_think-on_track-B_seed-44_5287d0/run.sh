#!/usr/bin/env bash
set -euo pipefail

# Ensure reference is indexed
if [ ! -f data/ref/chrM.fa.fai ]; then
    samtools faidx data/ref/chrM.fa
fi

# List of samples
samples=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Process each sample
for sample in "${samples[@]}"; do
    # Check if BAM and BAM.bai exist
    if [ ! -f results/${sample}.bam ] || [ ! -f results/${sample}.bam.bai ]; then
        # Align with bwa
        bwa mem -t 4 data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
        samtools view -bS - | \
        samtools sort -o results/${sample}.bam -@ 4
        samtools index results/${sample}.bam
    fi

    # Check if VCF.gz and TBI exist
    if [ ! -f results/${sample}.vcf.gz ] || [ ! -f results/${sample}.vcf.gz.tbi ]; then
        # Call variants with bcftools
        bcftools mpileup -f data/ref/chrM.fa results/${sample}.bam | \
        bcftools call -mv -a -o results/${sample}.vcf.gz
        tabix -p vcf results/${sample}.vcf.gz
    fi
done

# Generate TSV files for each sample
for sample in "${samples[@]}"; do
    if [ ! -f results/${sample}.tsv ]; then
        # Extract variants without header
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz | \
        grep -v '^#' | \
        awk -v sample="$sample" '{print sample "\t" $0}' > results/${sample}.tsv
    fi
done

# Combine all TSVs into collapsed.tsv
if [ ! -f results/collapsed.tsv ]; then
    # Combine all TSVs and add header
    cat results/*.tsv | sed '1i\sample\tchrom\tpos\tref\talt\taf' > results/collapsed.tsv
fi