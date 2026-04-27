#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF="data/ref/chrM.fa"
RAW_DIR="data/raw"
OUT_DIR="results"

mkdir -p "${OUT_DIR}"

# 2. Reference indexing (idempotent)
if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.bwt" ]]; then
    samtools faidx "${REF}"
    bwa index "${REF}"
fi

# 3-7. Per-sample alignment, sorting, indexing, calling, compression
for SAMPLE in "${SAMPLES[@]}"; do
    BAM="${OUT_DIR}/${SAMPLE}.bam"
    BAI="${OUT_DIR}/${SAMPLE}.bam.bai"
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    TBI="${OUT_DIR}/${SAMPLE}.vcf.gz.tbi"

    # Idempotency guard: if final VCF index exists, skip sample
    if [[ -f "${TBI}" ]]; then
        continue
    fi

    # 3 & 4. Alignment and sorting
    if [[ ! -f "${BAM}" ]]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
            "${REF}" \
            "${RAW_DIR}/${SAMPLE}_1.fq.gz" \
            "${RAW_DIR}/${SAMPLE}_2.fq.gz" | \
        samtools sort -@ "${THREADS}" -o "${BAM}" -
    fi

    # 5. BAM indexing
    if [[ ! -f "${BAI}" ]]; then
        samtools index -@ "${THREADS}" "${BAM}"
    fi

    # 6. Variant calling
    VCF_TMP="${OUT_DIR}/${SAMPLE}.vcf"
    if [[ ! -f "${VCF_GZ}" ]]; then
        lofreq call-parallel \
            --pp-threads "${THREADS}" \
            -f "${REF}" \
            -o "${VCF_TMP}" \
            "${BAM}"
    fi

    # 7. Compression and indexing
    if [[ ! -f "${VCF_GZ}" ]]; then
        bgzip -c "${VCF_TMP}" > "${VCF_GZ}"
        tabix -p vcf "${VCF_GZ}"
        rm -f "${VCF_TMP}"
    fi
done

# 8. Collapse step
COLLAPSED="${OUT_DIR}/collapsed.tsv"
NEED_COLLAPSE=false

# Check if collapsed.tsv exists
if [[ ! -f "${COLLAPSED}" ]]; then
    NEED_COLLAPSE=true
else
    # Check if any VCF is newer than collapsed.tsv
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
        if [[ "${VCF_GZ}" -nt "${COLLAPSED}" ]]; then
            NEED_COLLAPSE=true
            break
        fi
    done
fi

if [[ "${NEED_COLLAPSE}" == true ]]; then
    {
        printf "sample\tchrom\tpos\tref\talt\taf\n"
        for SAMPLE in "${SAMPLES[@]}"; do
            VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
            bcftools query -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "${VCF_GZ}"
        done
    } > "${COLLAPSED}"
fi