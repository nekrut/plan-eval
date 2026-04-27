#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF_FA="data/ref/chrM.fa"
RAW_DIR="data/raw"
RESULTS_DIR="results"

mkdir -p "${RESULTS_DIR}"

# 1. Reference indexing - BWA
if [[ ! -f "${REF_FA}.bwt" ]]; then
    bwa index "${REF_FA}"
fi

# 2. Reference indexing - samtools faidx
if [[ ! -f "${REF_FA}.fai" ]]; then
    samtools faidx "${REF_FA}"
fi

# 3. Per-sample alignment + sort
for sample in "${SAMPLES[@]}"; do
    BAM_OUT="${RESULTS_DIR}/${sample}.bam"
    if [[ ! -f "${BAM_OUT}" ]]; then
        echo "Processing alignment for ${sample}..."
        # Note: Using literal \t for tab character in RG header
        bwa mem -t "${THREADS}" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" "${REF_FA}" \
            "${RAW_DIR}/${sample}_1.fq.gz" "${RAW_DIR}/${sample}_2.fq.gz" | \
        samtools sort -@ "${THREADS}" -o "${BAM_OUT}" -
    fi
done

# 4. BAM index
for sample in "${SAMPLES[@]}"; do
    BAM_OUT="${RESULTS_DIR}/${sample}.bam"
    BAI_OUT="${RESULTS_DIR}/${sample}.bam.bai"
    if [[ ! -f "${BAI_OUT}" ]]; then
        samtools index -@ "${THREADS}" "${BAM_OUT}"
    fi
done

# 5. Variant calling - LoFreq
for sample in "${SAMPLES[@]}"; do
    VCF_OUT_BASE="${RESULTS_DIR}/${sample}"
    VCF_OUT_GZ="${VCF_OUT_BASE}.vcf.gz"
    
    # Check for existence of the final indexed file (sentinel)
    if [[ ! -f "${VCF_OUT_GZ}.tbi" ]]; then
        echo "Calling variants for ${sample}..."
        # lofreq call-parallel creates .vcf, which we then compress/index
        lofreq call-parallel --pp-threads "${THREADS}" -f "${REF_FA}" -o "${VCF_OUT_BASE}.vcf" "${RESULTS_DIR}/${sample}.bam"
    fi
done

# 6. VCF compression + tabix index
for sample in "${SAMPLES[@]}"; do
    VCF_OUT_BASE="${RESULTS_DIR}/${sample}"
    VCF_OUT_GZ="${VCF_OUT_BASE}.vcf.gz"
    
    if [[ ! -f "${VCF_OUT_GZ}.tbi" ]]; then
        echo "Compressing and indexing VCF for ${sample}..."
        # bgzip operates in place, deleting .vcf
        bgzip -f "${VCF_OUT_BASE}.vcf"
        tabix -p vcf "${VCF_OUT_GZ}"
    fi
done

# 7. Collapsed TSV (Always rebuild)
echo "Generating collapsed summary table..."
# Header
printf "sample\tchrom\tpos\tref\talt\taf\n" > "${RESULTS_DIR}/collapsed.tsv"

# Per sample append
for sample in "${SAMPLES[@]}"; do
    VCF_OUT_GZ="${RESULTS_DIR}/${sample}.vcf.gz"
    if [[ -f "${VCF_OUT_GZ}" ]]; then
        # bcftools query format: %CHROM\t%POS\t%REF\t%ALT\t%INFO/AF
        # awk prepends the sample name ($s) to the tab-separated output from bcftools
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${VCF_OUT_GZ}" | \
        awk -v s="${sample}" 'BEGIN{OFS="\t"}{print s,$0}' >> "${RESULTS_DIR}/collapsed.tsv"
    fi
done

exit 0