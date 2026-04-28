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

# Steps 3-7: Per-sample alignment, sorting, indexing, variant calling
for SAMPLE in ${SAMPLES}; do
    FQ1="${RAW}/${SAMPLE}_1.fq.gz"
    FQ2="${RAW}/${SAMPLE}_2.fq.gz"
    BAM="${RESULTS}/${SAMPLE}.bam"
    BAI="${RESULTS}/${SAMPLE}.bam.bai"
    VCF="${RESULTS}/${SAMPLE}.vcf"
    VCFGZ="${RESULTS}/${SAMPLE}.vcf.gz"
    TBI="${RESULTS}/${SAMPLE}.vcf.gz.tbi"

    # Step 3+4: Align and sort into BAM
    if [ ! -f "${BAM}" ]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
            "${REF}" "${FQ1}" "${FQ2}" \
            | samtools sort -@ "${THREADS}" -o "${BAM}"
    fi

    # Step 5: Index BAM
    if [ ! -f "${BAI}" ]; then
        samtools index -@ "${THREADS}" "${BAM}"
    fi

    # Steps 6+7: Call variants, compress, index
    if [ ! -f "${TBI}" ]; then
        if [ ! -f "${VCFGZ}" ]; then
            if [ ! -f "${VCF}" ]; then
                lofreq call-parallel --pp-threads "${THREADS}" --verbose \
                    --ref "${REF}" \
                    --out "${VCF}" \
                    "${BAM}"
            fi
            bgzip "${VCF}"
        fi
        tabix -p vcf "${VCFGZ}"
    fi
done

# Step 8: Build collapsed.tsv
COLLAPSED="${RESULTS}/collapsed.tsv"
NEED_COLLAPSE=0

if [ ! -f "${COLLAPSED}" ]; then
    NEED_COLLAPSE=1
else
    for SAMPLE in ${SAMPLES}; do
        if [ "${RESULTS}/${SAMPLE}.vcf.gz" -nt "${COLLAPSED}" ]; then
            NEED_COLLAPSE=1
            break
        fi
    done
fi

if [ "${NEED_COLLAPSE}" -eq 1 ]; then
    {
        printf 'sample\tchrom\tpos\tref\talt\taf\n'
        for SAMPLE in ${SAMPLES}; do
            bcftools query \
                -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
                "${RESULTS}/${SAMPLE}.vcf.gz"
        done
    } > "${COLLAPSED}"
fi