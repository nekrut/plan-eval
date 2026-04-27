#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF="data/ref/chrM.fa"

mkdir -p results

# 2. Reference indexing (once)
if [ ! -f "data/ref/chrM.fa.fai" ] || [ ! -f "data/ref/chrM.fa.bwt" ]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# 3-7. Per-sample processing
for sample in "${SAMPLES[@]}"; do
    bam="results/${sample}.bam"
    bai="results/${sample}.bam.bai"
    vcf_gz="results/${sample}.vcf.gz"
    vcf_tbi="results/${sample}.vcf.gz.tbi"
    vcf="results/${sample}.vcf"

    # 3-5. Alignment, sorting, indexing
    if [ ! -f "$bai" ] || [ "data/raw/${sample}_1.fq.gz" -nt "$bai" ] || [ "data/raw/${sample}_2.fq.gz" -nt "$bai" ]; then
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            "$REF" \
            "data/raw/${sample}_1.fq.gz" \
            "data/raw/${sample}_2.fq.gz" \
            | samtools sort -@ "$THREADS" -o "$bam" -
        samtools index -@ "$THREADS" "$bam"
    fi

    # 6-7. Variant calling, compression, indexing
    if [ ! -f "$vcf_tbi" ] || [ "$bam" -nt "$vcf_tbi" ]; then
        lofreq call-parallel -f "$REF" -i "$bam" -o "$vcf" --pp-threads "$THREADS"
        bgzip -c "$vcf" > "$vcf_gz"
        tabix -p vcf "$vcf_gz"
        rm -f "$vcf"
    fi
done

# 8. Collapse step
collapsed="results/collapsed.tsv"
need_rebuild=false
for sample in "${SAMPLES[@]}"; do
    if [ ! -f "results/${sample}.vcf.gz" ] || [ ! -f "results/${sample}.vcf.gz.tbi" ]; then
        need_rebuild=true
        break
    fi
    if [ ! -f "$collapsed" ] || [ "results/${sample}.vcf.gz" -nt "$collapsed" ]; then
        need_rebuild=true
        break
    fi
done

if $need_rebuild; then
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$collapsed"
    for sample in "${SAMPLES[@]}"; do
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "results/${sample}.vcf.gz" >> "$collapsed"
    done
fi