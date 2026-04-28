#!/usr/bin/env bash
set -euo pipefail

# Ensure results directory exists
mkdir -p results

# Index reference if not already done
if [ ! -f results/chrM.fa.bwt ]; then
    bwa index data/ref/chrM.fa
    samtools faidx data/ref/chrM.fa
    tabix -p fa data/ref/chrM.fa
fi

# Define samples
samples=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Loop through samples
for sample in "${samples[@]}"; do
    # Check if outputs already exist
    if [ -f "results/${sample}.vcf.gz" ]; then
        continue
    fi

    # 1. Map reads
    bwa mem -t 4 data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
    samtools view -b - > "results/${sample}.bam"

    # 2. Sort and index BAM
    samtools sort -@ 2 -o "results/${sample}.sorted.bam" "results/${sample}.bam"
    samtools index "results/${sample}.sorted.bam"
    mv "results/${sample}.sorted.bam" "results/${sample}.bam"

    # 3. Call variants with lofreq
    lofreq call-parallel -f data/ref/chrM.fa -o "results/${sample}.vcf" -m lofreq-mpileup -t 2 "results/${sample}.bam"

    # 4. Convert to bgzip and tabix
    bcftools view -Oz -o "results/${sample}.vcf.gz" "results/${sample}.vcf"
    tabix -p vcf "results/${sample}.vcf.gz"
done

# 5. Collapse VCFs
bcftools concat -Oz -o results/collapsed.vcf.gz results/*.vcf.gz
tabix -p vcf results/collapsed.vcf.gz

# 6. Extract and format columns
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' results/collapsed.vcf.gz | \
awk 'NR>1' > results/collapsed.tsv

# 7. Add header
echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv
cat results/collapsed.tsv > results/collapsed.tsv.tmp && mv results/collapsed.tsv.tmp results/collapsed.tsv

# 8. Cleanup intermediate files
rm -f results/*.vcf

exit 0