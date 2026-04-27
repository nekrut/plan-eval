#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF_FA="data/ref/chrM.fa"
RAW_DIR="data/raw"
RESULTS_DIR="results"

mkdir -p "$RESULTS_DIR"

# --- 1. Reference indexing — BWA ---
if [[ ! -f "${REF_FA}.bwt" ]]; then
    bwa index "$REF_FA"
fi

# --- 2. Reference indexing — samtools faidx ---
if [[ ! -f "${REF_FA}.fai" ]]; then
    samtools faidx "$REF_FA"
fi

# --- 3. Per-sample alignment + sort ---
for sample in "${SAMPLES[@]}"; do
    SAMPLE_BAM="${RESULTS_DIR}/${sample}.bam"
    if [[ ! -f "$SAMPLE_BAM" ]]; then
        # Construct the -R header string carefully to escape tabs
        RG_HEADER="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
        
        bwa mem -t "$THREADS" -R "$RG_HEADER" "$REF_FA" \
            "${RAW_DIR}/${sample}_1.fq.gz" "${RAW_DIR}/${sample}_2.fq.gz" | \
        samtools sort -@ "$THREADS" -o "$SAMPLE_BAM" -
    fi
done

# --- 4. BAM index ---
for sample in "${SAMPLES[@]}"; do
    SAMPLE_BAM="${RESULTS_DIR}/${sample}.bam"
    SAMPLE_BAI="${RESULTS_DIR}/${sample}.bam.bai"
    if [[ ! -f "$SAMPLE_BAI" ]]; then
        samtools index -@ "$THREADS" "$SAMPLE_BAM"
    fi
done

# --- 5. Variant calling — LoFreq ---
for sample in "${SAMPLES[@]}"; do
    SAMPLE_VCF="${RESULTS_DIR}/${sample}.vcf"
    SAMPLE_VCF_GZ="${RESULTS_DIR}/${sample}.vcf.gz"
    
    # Check if the final indexed file exists (Step 6 sentinel)
    if [[ ! -f "${RESULTS_DIR}/${sample}.vcf.gz.tbi" ]]; then
        lofreq call-parallel --pp-threads "$THREADS" \
            -f "$REF_FA" \
            -o "$SAMPLE_VCF" \
            "${RESULTS_DIR}/${sample}.bam"
    fi
done

# --- 6. VCF compression + tabix index ---
for sample in "${SAMPLES[@]}"; do
    SAMPLE_VCF="${RESULTS_DIR}/${sample}.vcf"
    SAMPLE_VCF_GZ="${RESULTS_DIR}/${sample}.vcf.gz"
    SAMPLE_VCF_TBI="${RESULTS_DIR}/${sample}.vcf.gz.tbi"

    if [[ ! -f "$SAMPLE_VCF_TBI" ]]; then
        # bgzip operates in place, creating .gz and deleting .vcf
        bgzip -f "$SAMPLE_VCF"
        tabix -p vcf "$SAMPLE_VCF.gz"
    fi
done

# --- 7. Collapsed TSV (Rebuild every run) ---
# Header
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${RESULTS_DIR}/collapsed.tsv"

# Per sample append
for sample in "${SAMPLES[@]}"; do
    SAMPLE_VCF_GZ="${RESULTS_DIR}/${sample}.vcf.gz"
    
    # bcftools query format: %CHROM\t%POS\t%REF\t%ALT\t%INFO/AF
    # awk prepends the sample name to the tab-separated output from bcftools
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$SAMPLE_VCF_GZ" | \
    awk -v s="$sample" 'BEGIN{OFS="\t"} {print s, $0}' >> "${RESULTS_DIR}/collapsed.tsv"
done