#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

mkdir -p results
cd "$(dirname "$0")"

# Reference indexing
if [[ ! -f data/ref/chrM.fa.amb ]]; then
    samtools faidx data/ref/chrM.fa
    bwa index data/ref/chrM.fa
fi

# Per-sample workflow
for sample in "${SAMPLES[@]}"; do
    bam="results/${sample}.bam"
    bai="results/${sample}.bam.bai"
    vcf="results/${sample}.vcf"
    vcf_gz="results/${sample}.vcf.gz"
    vcf_gz_tbi="results/${sample}.vcf.gz.tbi"

    # Skip if all outputs exist and are newer than inputs
    if [[ -f "$vcf_gz_tbi" ]] && [[ "$vcf_gz_tbi" -nt "data/raw/${sample}_1.fq.gz" ]] && [[ "$vcf_gz_tbi" -nt "data/raw/${sample}_2.fq.gz" ]]; then
        continue
    fi

    # Alignment
    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" data/ref/chrM.fa "data/raw/${sample}_1.fq.gz" "data/raw/${sample}_2.fq.gz" | \
    samtools sort -@ "$THREADS" -o "$bam"

    # Index BAM
    samtools index -@ "$THREADS" "$bam"

    # Variant calling
    lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o "$vcf" "$bam"

    # Compress and index VCF
    bgzip -c "$vcf" > "$vcf_gz"
    tabix -p vcf "$vcf_gz"
    rm -f "$vcf"
done

# Collapse step
if [[ ! -f results/collapsed.tsv ]] || [[ results/collapsed.tsv -lt results/M117-bl.vcf.gz.tbi ]] || [[ results/collapsed.tsv -lt results/M117-ch.vcf.gz.tbi ]] || [[ results/collapsed.tsv -lt results/M117C1-bl.vcf.gz.tbi ]] || [[ results/collapsed.tsv -lt results/M117C1-ch.vcf.gz.tbi ]]; then
    {
        echo -e "sample\tchrom\tpos\tref\talt\taf"
        for sample in "${SAMPLES[@]}"; do
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "results/${sample}.vcf.gz"
        done
    } > results/collapsed.tsv
fi