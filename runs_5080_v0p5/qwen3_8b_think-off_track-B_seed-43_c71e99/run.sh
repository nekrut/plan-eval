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
    lofreq call -f data/ref/chrM.fa -o results/${sample}.vcf.gz results/${sample}.bam || exit 1

    # Index VCF
    bcftools index -t z results/${sample}.vcf.gz || exit 1

    # Filter and collapse variants
    if [[ ! -f results/collapsed.tsv ]]; then
        awk 'NR==1 {print $1, $2, $3, $4, $5}' results/${sample}.vcf.gz | \
        awk '{count[$1"_"$2"_"$3"_"$4]++} END {for (key in count) {split(key, arr, "_"); print arr[1], arr[2], arr[3], arr[4], arr[5], count[key]}}' | \
        awk '{print $1, $2, $3, $4, $5, $6}' > results/collapsed.tsv || exit 1
    fi
done