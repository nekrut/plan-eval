#!/usr/bin/env bash
set -euo pipefail

mkdir -p results
find data/raw -type f -name "*.fq.gz" | while read -r fq1; do
    sample=$(basename "$fq1" | sed 's/_1.fq.gz//; s/_2.fq.gz//')
    fq2="${fq1/_1.fq.gz/_2.fq.gz}"
    if [[ ! -f "$fq2" ]]; then continue; fi

    bwa mem data/ref/chrM.fa "$fq1" "$fq2" | \
        samtools view -Sb - | \
        samtools sort -o results/"$sample".bam && \
        samtools index results/"$sample".bam

    lofreq call -f data/ref/chrM.fa -o results/"$sample".vcf.gz results/"$sample".bam
    bcftools index -t z results/"$sample".vcf.gz

    if [[ -f results/"$sample".vcf.gz ]]; then
        awk -v sample="$sample" '
            BEGIN { OFS="\t" }
            /^#/ { print; next }
            {
                split($1, a, "chr");
                chrom = a[2];
                pos = $4;
                ref = $5;
                alt = $7;
                af = $10;
                print sample, chrom, pos, ref, alt, af
            }
        ' results/"$sample".vcf.gz >> results/collapsed.tsv
    fi
done

if [[ -f results/collapsed.tsv ]]; then
    awk 'NR == 1 { print } NR > 1 && !seen[$1]++ { print }' results/collapsed.tsv > results/collapsed.tsv.tmp && mv results/collapsed.tsv.tmp results/collapsed.tsv
fi