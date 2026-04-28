#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

for sample in M117-bl M117-ch M117C1-bl M117C1-ch; do
    if [[ -f results/${sample}.bam && -f results/${sample}.bam.bai && -f results/${sample}.vcf.gz && -f results/${sample}.vcf.gz.tbi ]]; then
        continue
    fi

    # Align with bwa
    bwa mem data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools view -Sb - > results/${sample}.bam || exit 1

    # Index BAM
    samtools index results/${sample}.bam || exit 1

    # Call variants with lofreq
    lofreq call -f data/ref/chrM.fa -t 4 results/${sample}.bam | bcftools view -O z -o results/${sample}.vcf.gz || exit 1

    # Index VCF
    bcftools index results/${sample}.vcf.gz || exit 1

    # Extract variants and collapse
    if [[ ! -f results/collapsed.tsv ]]; then
        awk 'NR==1 {print $1, $2, $3, $4, $5}' results/${sample}.vcf.gz | cut -f1-5 | awk '{print $1, $2, $3, $4, $5}' > results/collapsed.tsv || exit 1
    else
        awk 'NR==1 {print $1, $2, $3, $4, $5}' results/${sample}.vcf.gz | cut -f1-5 | awk '{print $1, $2, $3, $4, $5}' >> results/collapsed.tsv || exit 1
    fi
done

# Add AF column
cut -f1-5 results/collapsed.tsv | awk '{print $1, $2, $3, $4, $5, 1}' > results/collapsed.tsv.tmp && mv results/collapsed.tsv.tmp results/collapsed.tsv || exit 1

# Sort and deduplicate
sort -k1,1 -k2,2n -k3,3n results/collapsed.tsv | uniq > results/collapsed.tsv.tmp && mv results/collapsed.tsv.tmp results/collapsed.tsv || exit 1

exit 0