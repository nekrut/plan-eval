#!/usr/bin/env bash
set -euo pipefail

mkdir -p results
find data/raw/ -name "*.fq.gz" | while read -r fq1; do
    fq2="${fq1}_2.fq.gz"
    sample=$(basename "$fq1" .fq.gz)
    if [[ ! -f "$fq2" ]]; then
        continue
    fi

    if [[ -f results/"$sample".bam ]]; then
        continue
    fi

    bwa mem data/ref/chrM.fa "$fq1" "$fq2" | samtools sort -@ 4 -o results/"$sample".bam
    samtools index results/"$sample".bam

    bcftools mpileup -d 'all' -f data/ref/chrM.fa -I -Ou results/"$sample".bam | \
        bcftools call -mv -O z -o results/"$sample".vcf.gz

    bcftools index -t results/"$sample".vcf.gz

    snpEff -v -c data/ref/chrM.fa -s results/"$sample".vcf.gz | \
        SnpSift filter "ANN[*].ALLELE[*].REF != '.' && ANN[*].ALLELE[*].ALT != '.' && ANN[*].ALLELE[*].REF != ANN[*].ALLELE[*].ALT" | \
        awk -F'\t' 'NR==1 {next} {print $1 "\t" $7 "\t" $10 "\t" $11 "\t" $12 "\t" $13}' | \
        cut -f1,2,3,4,5,6 | \
        sort -k1,1 -k2,2n -k3,3n | \
        uniq > results/collapsed.tsv
done

if [[ ! -f results/collapsed.tsv ]]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv
fi

sort -k1,1 -k2,2n -k3,3n results/collapsed.tsv | uniq > results/collapsed.tsv