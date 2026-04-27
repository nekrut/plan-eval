#!/usr/bin/env bash
set -euo pipefail
THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF=data/ref/chrM.fa
RESULTS=results
mkdir -p "$RESULTS"

# Reference indexing (idempotent)
if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.amb" ]]; then
    samtools faidx "$REF" -@ "$THREADS"
    bwa index "$REF" -p
fi

# Per-sample processing
for sample in "${SAMPLES[@]}"; do
    f1="data/raw/${sample}_1.fq.gz"
    f2="data/raw/${sample}_2.fq.gz"
    outbam="${RESULTS}/${sample}.bam"
    outbai="${RESULTS}/${sample}.bam.bai"
    outvcf="${RESULTS}/${sample}.vcf.gz"
    outtbi="${RESULTS}/${sample}.vcf.gz.tbi"

    # Alignment (idempotent: skip if bam exists)
    if [[ ! -f "$outbam" ]]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" "$REF" "$f1" "$f2" 2>/dev/null | \
        samtools sort -@ "$THREADS" -o "$outbam"
    fi

    # BAM indexing (idempotent)
    if [[ ! -f "$outbai" ]]; then
        samtools index -@ "$THREADS" "$outbam"
    fi

    # Variant calling (idempotent: skip if vcf.gz exists)
    if [[ ! -f "$outvcf" ]]; then
        lofreq call-parallel -f "$REF" -o "${RESULTS}/${sample}.vcf" --pp-threads "$THREADS" "$outbam" 2>/dev/null
        bgzip -f "${RESULTS}/${sample}.vcf"
        tabix -p vcf "$outvcf"
        rm -f "${RESULTS}/${sample}.vcf"
    fi
done

# Collapse step (idempotent: rebuild only if any input VCF is newer)
if [[ -f "${RESULTS}/collapsed.tsv" ]]; then
    newest_vcf=""
    for sample in "${SAMPLES[@]}"; do
        vcf="${RESULTS}/${sample}.vcf.gz"
        if [[ -f "$vcf" ]]; then
            if [[ -z "$newest_vcf" ]] || [[ "$vcf" -nt "$newest_vcf" ]]; then
                newest_vcf="$vcf"
            fi
        fi
    done
    if [[ -z "$newest_vcf" ]] || [[ "$newest_vcf" -nt "${RESULTS}/collapsed.tsv" ]]; then
        : > "${RESULTS}/collapsed.tsv"
        for sample in "${SAMPLES[@]}"; do
            bcftools query -f '%sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${RESULTS}/${sample}.vcf.gz" >> "${RESULTS}/collapsed.tsv"
        done
    fi
else
    : > "${RESULTS}/collapsed.tsv"
    for sample in "${SAMPLES[@]}"; do
        bcftools query -f '%sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${RESULTS}/${sample}.vcf.gz" >> "${RESULTS}/collapsed.tsv"
    done
fi