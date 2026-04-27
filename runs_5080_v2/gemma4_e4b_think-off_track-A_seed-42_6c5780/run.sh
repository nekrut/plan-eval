#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
mkdir -p results

# --- 1. Reference indexing — BWA ---
REF_FA="data/ref/chrM.fa"
if [[ ! -f "${REF_FA}.bwt" ]]; then
    bwa index "${REF_FA}"
fi

# --- 2. Reference indexing — samtools faidx ---
if [[ ! -f "${REF_FA}.fai" ]]; then
    samtools faidx "${REF_FA}"
fi

# --- 3. Per-sample alignment + sort ---
for sample in "${SAMPLES[@]}"; do
    BAM_OUT="results/${sample}.bam"
    if [[ ! -f "${BAM_OUT}" ]]; then
        # Constructing the RG header string carefully to escape tabs
        RG_HEADER="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
        
        bwa mem -t ${THREADS} -R "${RG_HEADER}" "${REF_FA}" \
            "data/raw/${sample}_1.fq.gz" "data/raw/${sample}_2.fq.gz" | \
        samtools sort -@ ${THREADS} -o "${BAM_OUT}" -
    fi
done

# --- 4. BAM index ---
for sample in "${SAMPLES[@]}"; do
    BAM_OUT="results/${sample}.bam"
    BAI_OUT="results/${sample}.bam.bai"
    if [[ ! -f "${BAI_OUT}" ]]; then
        samtools index -@ ${THREADS} "${BAM_OUT}"
    fi
done

# --- 5. Variant calling — LoFreq ---
for sample in "${SAMPLES[@]}"; do
    VCF_OUT="results/${sample}.vcf"
    VCF_GZ_OUT="results/${sample}.vcf.gz"
    
    # Check for existence of the final indexed file (.tbi) to determine if work is needed
    if [[ ! -f "${VCF_GZ_OUT}.tbi" ]]; then
        lofreq call-parallel --pp-threads ${THREADS} \
            -f "${REF_FA}" \
            -o "${VCF_OUT}" \
            results/${sample}.bam
    fi
done

# --- 6. VCF compression + tabix index ---
for sample in "${SAMPLES[@]}"; do
    VCF_OUT="results/${sample}.vcf"
    VCF_GZ_OUT="results/${sample}.vcf.gz"
    
    if [[ ! -f "${VCF_GZ_OUT}.tbi" ]]; then
        # bgzip operates in place, creating .gz and deleting uncompressed .vcf
        bgzip -f "${VCF_OUT}"
        tabix -p vcf "${VCF_GZ_OUT}"
    fi
done

# --- 7. Collapsed TSV (rebuild every run) ---
COLLAPSED_TSV="results/collapsed.tsv"

# Write header
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${COLLAPSED_TSV}"

# Append data for each sample
for sample in "${SAMPLES[@]}"; do
    VCF_GZ_OUT="results/${sample}.vcf.gz"
    
    # bcftools query format: %CHROM\t%POS\t%REF\t%ALT\t%INFO/AF
    # awk prepends the sample name, ensuring the output is tab-separated
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${VCF_GZ_OUT}" | \
    awk -v s="${sample}" 'BEGIN{OFS="\t"}{print s,$0}' >> "${COLLAPSED_TSV}"
done

exit 0