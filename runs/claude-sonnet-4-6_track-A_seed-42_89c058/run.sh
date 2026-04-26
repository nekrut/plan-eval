#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF="data/ref/chrM.fa"
RAW="data/raw"
RESULTS="results"

mkdir -p "${RESULTS}"

# Step 2: Reference indexing
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "${REF}"
fi

if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "${REF}"
fi

for SAMPLE in "${SAMPLES[@]}"; do
    BAM="${RESULTS}/${SAMPLE}.bam"
    BAI="${RESULTS}/${SAMPLE}.bam.bai"
    VCF_GZ="${RESULTS}/${SAMPLE}.vcf.gz"
    VCF_TBI="${RESULTS}/${SAMPLE}.vcf.gz.tbi"
    VCF_TMP="${RESULTS}/${SAMPLE}.vcf"
    FQ1="${RAW}/${SAMPLE}_1.fq.gz"
    FQ2="${RAW}/${SAMPLE}_2.fq.gz"

    # Step 3 & 4: Alignment and sort
    if [[ ! -f "${BAM}" ]]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
            "${REF}" "${FQ1}" "${FQ2}" \
            | samtools sort -@ "${THREADS}" -o "${BAM}"
    fi

    # Step 5: BAM indexing
    if [[ ! -f "${BAI}" ]]; then
        samtools index -@ "${THREADS}" "${BAM}"
    fi

    # Step 6 & 7: Variant calling, compression, indexing
    if [[ ! -f "${VCF_TBI}" ]]; then
        lofreq call-parallel --pp-threads "${THREADS}" \
            -f "${REF}" \
            -o "${VCF_TMP}" \
            "${BAM}"
        bgzip -f "${VCF_TMP}"
        tabix -p vcf "${VCF_GZ}"
        [[ -f "${VCF_TMP}" ]] && rm -f "${VCF_TMP}"
    fi
done

# Step 8: Collapse step
REBUILD_COLLAPSED=0
COLLAPSED="${RESULTS}/collapsed.tsv"

if [[ ! -f "${COLLAPSED}" ]]; then
    REBUILD_COLLAPSED=1
else
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF_GZ="${RESULTS}/${SAMPLE}.vcf.gz"
        if [[ "${VCF_GZ}" -nt "${COLLAPSED}" ]]; then
            REBUILD_COLLAPSED=1
            break
        fi
    done
fi

if [[ "${REBUILD_COLLAPSED}" -eq 1 ]]; then
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${COLLAPSED}"
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF_GZ="${RESULTS}/${SAMPLE}.vcf.gz"
        bcftools query \
            -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
            "${VCF_GZ}" >> "${COLLAPSED}"
    done
fi