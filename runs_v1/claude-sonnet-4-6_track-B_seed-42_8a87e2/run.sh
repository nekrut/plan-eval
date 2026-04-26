#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF="data/ref/chrM.fa"
RAW="data/raw"
RESULTS="results"

SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

mkdir -p "${RESULTS}"

# Index reference for bwa
if [ ! -f "${REF}.bwt" ]; then
    bwa index "${REF}"
fi

# Index reference for samtools / lofreq
if [ ! -f "${REF}.fai" ]; then
    samtools faidx "${REF}"
fi

for SAMPLE in "${SAMPLES[@]}"; do
    FQ1="${RAW}/${SAMPLE}_1.fq.gz"
    FQ2="${RAW}/${SAMPLE}_2.fq.gz"
    BAM="${RESULTS}/${SAMPLE}.bam"
    IQ_BAM="${RESULTS}/${SAMPLE}.iq.bam"
    VCF_GZ="${RESULTS}/${SAMPLE}.vcf.gz"

    # 1. Align and sort
    if [ ! -f "${BAM}" ]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
            "${REF}" "${FQ1}" "${FQ2}" \
        | samtools sort -@ "${THREADS}" -o "${BAM}"
    fi

    # 2. Index final BAM (required output)
    if [ ! -f "${BAM}.bai" ]; then
        samtools index "${BAM}"
    fi

    # 3. Add indel quality scores required by lofreq --call-indels
    if [ ! -f "${IQ_BAM}" ]; then
        lofreq indelqual --dindel -f "${REF}" "${BAM}" -o "${IQ_BAM}"
    fi

    if [ ! -f "${IQ_BAM}.bai" ]; then
        samtools index "${IQ_BAM}"
    fi

    # 4. Call variants (SNVs + indels) with lofreq
    if [ ! -f "${VCF_GZ}" ]; then
        TMPVCF="${RESULTS}/${SAMPLE}.tmp.vcf"
        lofreq call \
            -f "${REF}" \
            --call-indels \
            -o "${TMPVCF}" \
            "${IQ_BAM}"
        bcftools view -O z -o "${VCF_GZ}" "${TMPVCF}"
        rm -f "${TMPVCF}"
    fi

    # 5. Index VCF (required output)
    if [ ! -f "${VCF_GZ}.tbi" ]; then
        tabix -p vcf "${VCF_GZ}"
    fi
done

# Collapse all per-sample VCFs into one table
COLLAPSED="${RESULTS}/collapsed.tsv"
if [ ! -f "${COLLAPSED}" ]; then
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${COLLAPSED}"
    for SAMPLE in "${SAMPLES[@]}"; do
        bcftools query \
            -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
            "${RESULTS}/${SAMPLE}.vcf.gz" \
        >> "${COLLAPSED}"
    done
fi