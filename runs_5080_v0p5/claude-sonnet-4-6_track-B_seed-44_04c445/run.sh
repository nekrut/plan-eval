#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF="data/ref/chrM.fa"
RAW="data/raw"
RESULTS="results"
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

mkdir -p "${RESULTS}"

# ── Reference indexing ────────────────────────────────────────────────────────
[[ -f "${REF}.bwt" ]] || bwa index "${REF}"
[[ -f "${REF}.fai" ]] || samtools faidx "${REF}"

# ── Per-sample pipeline ───────────────────────────────────────────────────────
for SAMPLE in "${SAMPLES[@]}"; do
    BAM="${RESULTS}/${SAMPLE}.bam"
    IQ_BAM="${RESULTS}/${SAMPLE}.iq.bam"
    VCF_RAW="${RESULTS}/${SAMPLE}.raw.vcf"
    VCF_GZ="${RESULTS}/${SAMPLE}.vcf.gz"

    # 1. Align → sort → index
    if [[ ! -f "${BAM}.bai" ]]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
            "${REF}" \
            "${RAW}/${SAMPLE}_1.fq.gz" \
            "${RAW}/${SAMPLE}_2.fq.gz" \
        | samtools sort -@ "${THREADS}" -o "${BAM}"
        samtools index "${BAM}"
    fi

    # 2. Add indel qualities (required for lofreq indel calling)
    if [[ ! -f "${IQ_BAM}.bai" ]]; then
        lofreq indelqual --dindel -f "${REF}" -o "${IQ_BAM}" "${BAM}"
        samtools index "${IQ_BAM}"
    fi

    # 3. Variant calling → bgzip → tabix
    if [[ ! -f "${VCF_GZ}.tbi" ]]; then
        lofreq call-parallel --pp-threads "${THREADS}" \
            --call-indels \
            -f "${REF}" \
            -o "${VCF_RAW}" \
            "${IQ_BAM}"
        bcftools view -Oz -o "${VCF_GZ}" "${VCF_RAW}"
        tabix -p vcf "${VCF_GZ}"
        rm -f "${VCF_RAW}"
    fi
done

# ── Collapse all VCFs into a single table ─────────────────────────────────────
COLLAPSED="${RESULTS}/collapsed.tsv"
if [[ ! -f "${COLLAPSED}" ]]; then
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${COLLAPSED}"
    for SAMPLE in "${SAMPLES[@]}"; do
        bcftools query \
            -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
            "${RESULTS}/${SAMPLE}.vcf.gz" \
            >> "${COLLAPSED}"
    done
fi