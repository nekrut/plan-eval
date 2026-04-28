#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES="M117-bl M117-ch M117C1-bl M117C1-ch"
REF="data/ref/chrM.fa"
RAW="data/raw"
RESULTS="results"

mkdir -p "${RESULTS}"

# Step 2: Reference indexing
if [ ! -f "${REF}.fai" ]; then
    samtools faidx "${REF}"
fi

if [ ! -f "${REF}.bwt" ]; then
    bwa index "${REF}"
fi

# Steps 3-7: Per-sample alignment, sorting, indexing, variant calling, compression
for SAMPLE in ${SAMPLES}; do
    BAM="${RESULTS}/${SAMPLE}.bam"
    BAI="${RESULTS}/${SAMPLE}.bam.bai"
    VCF="${RESULTS}/${SAMPLE}.vcf"
    VCFGZ="${RESULTS}/${SAMPLE}.vcf.gz"
    TBI="${RESULTS}/${SAMPLE}.vcf.gz.tbi"
    FQ1="${RAW}/${SAMPLE}_1.fq.gz"
    FQ2="${RAW}/${SAMPLE}_2.fq.gz"

    # Step 3+4: Alignment and sorting
    if [ ! -f "${BAM}" ]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
            "${REF}" "${FQ1}" "${FQ2}" \
            | samtools sort -@ "${THREADS}" -o "${BAM}"
    fi

    # Step 5: BAM indexing
    if [ ! -f "${BAI}" ]; then
        samtools index -@ "${THREADS}" "${BAM}"
    fi

    # Step 6: Variant calling
    if [ ! -f "${VCFGZ}" ] && [ ! -f "${VCF}" ]; then
        lofreq call-parallel --pp-threads "${THREADS}" --verbose \
            --ref "${REF}" \
            --out "${VCF}" \
            "${BAM}"
    fi

    # Step 7: Compression and indexing
    if [ ! -f "${TBI}" ]; then
        if [ ! -f "${VCFGZ}" ] && [ -f "${VCF}" ]; then
            bgzip "${VCF}"
        fi
        tabix -p vcf "${VCFGZ}"
        if [ -f "${VCF}" ]; then
            rm -f "${VCF}"
        fi
    fi
done

# Step 8: Collapse all samples into results/collapsed.tsv
NEED_COLLAPSE=0
if [ ! -f "${RESULTS}/collapsed.tsv" ]; then
    NEED_COLLAPSE=1
else
    for SAMPLE in ${SAMPLES}; do
        VCFGZ="${RESULTS}/${SAMPLE}.vcf.gz"
        if [ "${VCFGZ}" -nt "${RESULTS}/collapsed.tsv" ]; then
            NEED_COLLAPSE=1
            break
        fi
    done
fi

if [ "${NEED_COLLAPSE}" -eq 1 ]; then
    TMPTSV="${RESULTS}/collapsed.tsv.tmp"
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${TMPTSV}"
    for SAMPLE in ${SAMPLES}; do
        VCFGZ="${RESULTS}/${SAMPLE}.vcf.gz"
        bcftools query \
            -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
            "${VCFGZ}" >> "${TMPTSV}"
    done
    mv "${TMPTSV}" "${RESULTS}/collapsed.tsv"
fi