#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# 1. Reference indexing — BWA
if [ ! -f data/ref/chrM.fa.bwt ]; then
    bwa index data/ref/chrM.fa
fi

# 2. Reference indexing — samtools faidx
if [ ! -f data/ref/chrM.fa.fai ]; then
    samtools faidx data/ref/chrM.fa
fi

# Per-sample processing loop
for sample in "${SAMPLES[@]}"; do
    echo "Processing sample: ${sample}"

    # 3. Per-sample alignment + sort
    INPUT_R1="data/raw/${sample}_1.fq.gz"
    INPUT_R2="data/raw/${sample}_2.fq.gz"
    OUTPUT_BAM="results/${sample}.bam"

    if [ ! -f "${OUTPUT_BAM}" ]; then
        bwa mem -t ${THREADS} -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" data/ref/chrM.fa "${INPUT_R1}" "${INPUT_R2}" | \
        samtools sort -@ ${THREADS} -o "${OUTPUT_BAM}" -
    else
        : # Idempotent: Skip if BAM exists
    fi

    # 4. BAM index
    OUTPUT_BAI="${OUTPUT_BAM}.bai"
    if [ ! -f "${OUTPUT_BAI}" ]; then
        samtools index -@ ${THREADS} "${OUTPUT_BAM}"
    fi

    # 5. Variant calling — LoFreq
    OUTPUT_VCF="results/${sample}.vcf"
    if [ ! -f "${OUTPUT_VCF}" ]; then
        lofreq call-parallel --pp-threads ${THREADS} -f data/ref/chrM.fa -o "${OUTPUT_VCF}" "${OUTPUT_BAM}"
    fi

    # 6. VCF compression + tabix index
    OUTPUT_VCF_GZ="results/${sample}.vcf.gz"
    if [ ! -f "${OUTPUT_VCF_GZ}" ]; then
        bgzip -f "${OUTPUT_VCF}"
        mv "${OUTPUT_VCF}" "${OUTPUT_VCF_GZ}.tmp" # Rename to avoid issues if bgzip fails partially
        mv "${OUTPUT_VCF_GZ}.tmp" "${OUTPUT_VCF_GZ}"
    fi

    if [ ! -f "${OUTPUT_VCF_GZ}.tbi" ]; then
        tabix -p vcf "${OUTPUT_VCF_GZ}"
    fi
done

# 7. Collapsed TSV
COLLAPSED_TSV="results/collapsed.tsv"

# Initialize header
if [ ! -f "${COLLAPSED_TSV}" ]; then
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${COLLAPSED_TSV}"
fi

# Append data for each sample
for sample in "${SAMPLES[@]}"; do
    OUTPUT_VCF_GZ="results/${sample}.vcf.gz"
    if [ -f "${OUTPUT_VCF_GZ}" ]; then
        # Use bcftools query to extract fields and awk to prepend sample name
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${OUTPUT_VCF_GZ}" | \
        awk -v s="${sample}" 'BEGIN{OFS="\t"} {print s, $0}' >> "${COLLAPSED_TSV}"
    fi
done

exit 0