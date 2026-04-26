#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES="M117-bl M117-ch M117C1-bl M117C1-ch"
REF="data/ref/chrM.fa"

mkdir -p results

# Step 2: Reference indexing
if [ ! -f "${REF}.fai" ]; then
    samtools faidx "${REF}"
fi

if [ ! -f "data/ref/chrM.fa.bwt" ]; then
    bwa index "${REF}"
fi

# Steps 3-7: Per-sample processing
for SAMPLE in ${SAMPLES}; do
    FQ1="data/raw/${SAMPLE}_1.fq.gz"
    FQ2="data/raw/${SAMPLE}_2.fq.gz"
    BAM="results/${SAMPLE}.bam"
    BAI="results/${SAMPLE}.bam.bai"
    VCF_TMP="results/${SAMPLE}.vcf"
    VCF_GZ="results/${SAMPLE}.vcf.gz"
    VCF_TBI="results/${SAMPLE}.vcf.gz.tbi"

    # Step 3+4: Align and sort
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

    # Step 6: Variant calling
    if [ ! -f "${VCF_GZ}" ]; then
        lofreq call-parallel --pp-threads "${THREADS}" \
            -f "${REF}" \
            -o "${VCF_TMP}" \
            "${BAM}"

        # Step 7: Compress and index
        bgzip -c "${VCF_TMP}" > "${VCF_GZ}"
        tabix -p vcf "${VCF_GZ}"
        rm -f "${VCF_TMP}"
    fi

    # Ensure index exists even if VCF_GZ was already present
    if [ ! -f "${VCF_TBI}" ]; then
        tabix -p vcf "${VCF_GZ}"
    fi
done

# Step 8: Collapse
REBUILD_TSV=0
for SAMPLE in ${SAMPLES}; do
    VCF_GZ="results/${SAMPLE}.vcf.gz"
    if [ ! -f "results/collapsed.tsv" ] || [ "${VCF_GZ}" -nt "results/collapsed.tsv" ]; then
        REBUILD_TSV=1
        break
    fi
done

if [ ! -f "results/collapsed.tsv" ]; then
    REBUILD_TSV=1
fi

if [ "${REBUILD_TSV}" -eq 1 ]; then
    TMPTSV="results/collapsed.tsv.tmp"
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${TMPTSV}"
    for SAMPLE in ${SAMPLES}; do
        VCF_GZ="results/${SAMPLE}.vcf.gz"
        bcftools query \
            -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
            "${VCF_GZ}" >> "${TMPTSV}"
    done
    mv "${TMPTSV}" "results/collapsed.tsv"
fi