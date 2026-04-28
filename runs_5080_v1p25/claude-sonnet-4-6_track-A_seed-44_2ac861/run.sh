#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES="M117-bl M117-ch M117C1-bl M117C1-ch"
REF="data/ref/chrM.fa"
RAW="data/raw"
RESULTS="results"

mkdir -p "${RESULTS}"

# Step 2: Reference indexing
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "${REF}"
fi

if [[ ! -f "data/ref/chrM.fa.bwt" ]]; then
    bwa index "${REF}"
fi

# Steps 3-7: Per-sample processing
for SAMPLE in ${SAMPLES}; do
    FQ1="${RAW}/${SAMPLE}_1.fq.gz"
    FQ2="${RAW}/${SAMPLE}_2.fq.gz"
    BAM="${RESULTS}/${SAMPLE}.bam"
    BAI="${RESULTS}/${SAMPLE}.bam.bai"
    VCF="${RESULTS}/${SAMPLE}.vcf"
    VCFGZ="${RESULTS}/${SAMPLE}.vcf.gz"
    TBI="${RESULTS}/${SAMPLE}.vcf.gz.tbi"

    # Steps 3-4: Alignment → sorted BAM
    if [[ ! -f "${BAM}" ]]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
            "${REF}" "${FQ1}" "${FQ2}" | \
            samtools sort -@ "${THREADS}" -o "${BAM}"
    fi

    # Step 5: BAM indexing
    if [[ ! -f "${BAI}" ]]; then
        samtools index -@ "${THREADS}" "${BAM}"
    fi

    # Step 6: Variant calling
    if [[ ! -f "${VCF}" ]] && [[ ! -f "${VCFGZ}" ]]; then
        lofreq call-parallel --pp-threads "${THREADS}" \
            -f "${REF}" \
            -o "${VCF}" \
            "${BAM}"
    fi

    # Step 7: VCF compression and indexing
    if [[ ! -f "${TBI}" ]]; then
        if [[ ! -f "${VCFGZ}" ]]; then
            bgzip "${VCF}"
        fi
        tabix -p vcf "${VCFGZ}"
        if [[ -f "${VCF}" ]]; then
            rm "${VCF}"
        fi
    fi
done

# Step 8: Collapse
COLLAPSED="${RESULTS}/collapsed.tsv"
REBUILD=0

if [[ ! -f "${COLLAPSED}" ]]; then
    REBUILD=1
else
    for SAMPLE in ${SAMPLES}; do
        if [[ "${RESULTS}/${SAMPLE}.vcf.gz" -nt "${COLLAPSED}" ]]; then
            REBUILD=1
            break
        fi
    done
fi

if [[ "${REBUILD}" -eq 1 ]]; then
    {
        printf 'sample\tchrom\tpos\tref\talt\taf\n'
        for SAMPLE in ${SAMPLES}; do
            bcftools query \
                -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
                "${RESULTS}/${SAMPLE}.vcf.gz"
        done
    } > "${COLLAPSED}"
fi