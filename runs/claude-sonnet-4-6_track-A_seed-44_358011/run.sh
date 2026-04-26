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

if [ ! -f "${REF}.bwt" ]; then
    bwa index "${REF}"
fi

# Steps 3-7: Per-sample alignment, sorting, indexing, variant calling
for SAMPLE in ${SAMPLES}; do
    BAM="results/${SAMPLE}.bam"
    BAI="results/${SAMPLE}.bam.bai"
    VCF_GZ="results/${SAMPLE}.vcf.gz"
    VCF_TBI="results/${SAMPLE}.vcf.gz.tbi"
    VCF_TMP="results/${SAMPLE}.vcf"
    R1="data/raw/${SAMPLE}_1.fq.gz"
    R2="data/raw/${SAMPLE}_2.fq.gz"

    # Steps 3-4: Align and sort
    if [ ! -f "${BAM}" ]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
            "${REF}" "${R1}" "${R2}" \
            | samtools sort -@ "${THREADS}" -o "${BAM}"
    fi

    # Step 5: Index BAM
    if [ ! -f "${BAI}" ]; then
        samtools index -@ "${THREADS}" "${BAM}"
    fi

    # Steps 6-7: Call variants, compress, index VCF
    if [ ! -f "${VCF_TBI}" ]; then
        if [ ! -f "${VCF_GZ}" ]; then
            lofreq call-parallel --pp-threads "${THREADS}" \
                -f "${REF}" \
                -o "${VCF_TMP}" \
                "${BAM}"
            bgzip -c "${VCF_TMP}" > "${VCF_GZ}"
            rm -f "${VCF_TMP}"
        fi
        tabix -p vcf "${VCF_GZ}"
    fi
done

# Step 8: Build collapsed.tsv
rebuild_collapsed() {
    {
        printf 'sample\tchrom\tpos\tref\talt\taf\n'
        for SAMPLE in ${SAMPLES}; do
            bcftools query \
                -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
                "results/${SAMPLE}.vcf.gz"
        done
    } > results/collapsed.tsv
}

if [ ! -f "results/collapsed.tsv" ]; then
    rebuild_collapsed
else
    needs_rebuild=0
    for SAMPLE in ${SAMPLES}; do
        if [ "results/${SAMPLE}.vcf.gz" -nt "results/collapsed.tsv" ]; then
            needs_rebuild=1
            break
        fi
    done
    if [ "${needs_rebuild}" -eq 1 ]; then
        rebuild_collapsed
    fi
fi