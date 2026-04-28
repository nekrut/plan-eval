#!/usr/bin/env bash
set -euo pipefail
THREADS=4
SAMPLES="M117-bl M117-ch M117C1-bl M117C1-ch"
REF="data/ref/chrM.fa"
RAW_PREFIX="data/raw"
RESULTS="results"

mkdir -p "$RESULTS"

REF_INDEXED=0
if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.amb" ]] || [[ ! -f "${REF}.ann" ]] || [[ ! -f "${REF}.bwt" ]] || [[ ! -f "${REF}.pac" ]] || [[ ! -f "${REF}.sa" ]]; then
    samtools faidx "$REF"
    bwa index "$REF"
    REF_INDEXED=1
fi

for SAMPLE in $SAMPLES; do
    SAMPLE_1="${RAW_PREFIX}/${SAMPLE}_1.fq.gz"
    SAMPLE_2="${RAW_PREFIX}/${SAMPLE}_2.fq.gz"
    BAM="${RESULTS}/${SAMPLE}.bam"
    BAM_BAI="${RESULTS}/${SAMPLE}.bam.bai"
    VCF="${RESULTS}/${SAMPLE}.vcf.gz"
    VCF_TBI="${RESULTS}/${SAMPLE}.vcf.gz.tbi"

    if [[ -f "$BAM" ]] && [[ -f "$BAM_BAI" ]] && [[ -f "$VCF" ]] && [[ -f "$VCF_TBI" ]]; then
        continue
    fi

    if [[ -f "$BAM" ]]; then
        rm -f "$BAM" "$BAM_BAI"
    fi

    bwa mem -t "$THREADS" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" "$REF" "$SAMPLE_1" "$SAMPLE_2" | samtools sort -@ "$THREADS" -o "$BAM"

    samtools index -@ "$THREADS" "$BAM"

    lofreq call-parallel --pp-threads "$THREADS" --verbose --ref "$REF" --out "$RESULTS/${SAMPLE}.vcf" "$BAM"

    bgzip -f "$RESULTS/${SAMPLE}.vcf"
    tabix -p vcf "$VCF"
    rm -f "$RESULTS/${SAMPLE}.vcf"
done

if [[ -f "${RESULTS}/collapsed.tsv" ]]; then
    for SAMPLE in $SAMPLES; do
        if [[ ! -f "${RESULTS}/${SAMPLE}.vcf.gz" ]]; then
            rm -f "${RESULTS}/collapsed.tsv"
            break
        fi
    done
fi

if [[ ! -f "${RESULTS}/collapsed.tsv" ]]; then
    {
        printf 'sample\tchrom\tpos\tref\talt\taf\n'
        for SAMPLE in $SAMPLES; do
            bcftools query -f '%sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${RESULTS}/${SAMPLE}.vcf.gz"
        done
    } > "${RESULTS}/collapsed.tsv"
fi