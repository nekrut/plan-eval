#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES="M117-bl M117-ch M117C1-bl M117C1-ch"
REF="data/ref/chrM.fa"

mkdir -p results

# Step 2: Reference indexing
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "${REF}"
fi

if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "${REF}"
fi

# Steps 3-7: Per-sample processing
for SAMPLE in ${SAMPLES}; do
    R1="data/raw/${SAMPLE}_1.fq.gz"
    R2="data/raw/${SAMPLE}_2.fq.gz"
    BAM="results/${SAMPLE}.bam"
    BAI="results/${SAMPLE}.bam.bai"
    VCF="results/${SAMPLE}.vcf"
    VCFGZ="results/${SAMPLE}.vcf.gz"
    TBI="results/${SAMPLE}.vcf.gz.tbi"

    # Step 3-4: Align and sort
    if [[ ! -f "${BAM}" ]]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
            "${REF}" "${R1}" "${R2}" \
            | samtools sort -@ "${THREADS}" -o "${BAM}"
    fi

    # Step 5: Index BAM
    if [[ ! -f "${BAI}" ]]; then
        samtools index -@ "${THREADS}" "${BAM}"
    fi

    # Step 6: Variant calling
    if [[ ! -f "${VCFGZ}" ]]; then
        lofreq call-parallel --pp-threads "${THREADS}" \
            -f "${REF}" \
            -o "${VCF}" \
            "${BAM}"

        # Step 7: Compress and index
        bgzip -@ "${THREADS}" "${VCF}"
        tabix -p vcf "${VCFGZ}"
    elif [[ ! -f "${TBI}" ]]; then
        tabix -p vcf "${VCFGZ}"
    fi
done

# Step 8: Collapse all samples into results/collapsed.tsv
NEWEST_VCF=""
for SAMPLE in ${SAMPLES}; do
    VCFGZ="results/${SAMPLE}.vcf.gz"
    if [[ -z "${NEWEST_VCF}" ]] || [[ "${VCFGZ}" -nt "${NEWEST_VCF}" ]]; then
        NEWEST_VCF="${VCFGZ}"
    fi
done

if [[ ! -f "results/collapsed.tsv" ]] || [[ "${NEWEST_VCF}" -nt "results/collapsed.tsv" ]]; then
    {
        printf 'sample\tchrom\tpos\tref\talt\taf\n'
        for SAMPLE in ${SAMPLES}; do
            bcftools query \
                -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
                "results/${SAMPLE}.vcf.gz"
        done
    } > results/collapsed.tsv
fi