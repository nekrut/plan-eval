#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF_FA="data/ref/chrM.fa"
RESULTS_DIR="results"
COLLAPSED_FILE="${RESULTS_DIR}/collapsed.tsv"

# --- 1. Setup and Idempotency Check ---

mkdir -p "${RESULTS_DIR}"

# Check if all required outputs exist to determine if work is needed
ALL_ARTIFACTS_EXIST=true
for sample in "${SAMPLES[@]}"; do
    if [[ ! -f "${RESULTS_DIR}/${sample}.vcf.gz.tbi" ]]; then
        ALL_ARTIFACTS_EXIST=false
        break
    fi
done

if $ALL_ARTIFACTS_EXIST && [[ -f "${COLLAPSED_FILE}" ]]; then
    # Basic check: If the final collapsed file exists, assume success for this run.
    # A more rigorous check would compare timestamps, but for simplicity and idempotency,
    # checking the final output is sufficient.
    exit 0
fi

# --- 2. Reference Indexing ---

echo "Indexing reference genome..." >&2
if [[ ! -f "${REF_FA}.fai" ]]; then
    samtools faidx "${REF_FA}"
fi

if [[ ! -f "${REF_FA}.bwt" ]]; then
    bwa index "${REF_FA}"
fi

# --- 3. Per-sample Alignment (bwa mem) ---

for sample in "${SAMPLES[@]}"; do
    echo "Processing sample: ${sample}" >&2
    R1="${data/raw/${sample}_1.fq.gz"
    R2="${data/raw/${sample}_2.fq.gz"
    
    # Check if input files exist before proceeding
    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        echo "Skipping ${sample}: Input FASTQ files not found." >&2
        continue
    fi

    # Read Group definition: -R "@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA"
    # Note: The literal backslash-t must be passed carefully.
    RG_LINE="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"

    # Check if BAM already exists and is up-to-date (simple check)
    if [[ -f "${RESULTS_DIR}/${sample}.bam" ]]; then
        # A full timestamp check is complex; we rely on the subsequent steps failing if inputs are missing.
        : # Assume we proceed unless the BAM is missing entirely.
    else
        echo "  -> Running bwa mem for ${sample}..." >&2
        bwa mem -t ${THREADS} -R "${RG_LINE}" "${REF_FA}" "${R1}" "${R2}" | \
        samtools sort -@ ${THREADS} -o "${RESULTS_DIR}/${sample}.bam"
    fi
done

# --- 4. BAM Indexing ---

for sample in "${SAMPLES[@]}"; do
    BAM_FILE="${RESULTS_DIR}/${sample}.bam"
    if [[ -f "${BAM_FILE}" && ! -f "${BAM_FILE}.bai" ]]; then
        echo "Indexing ${sample} BAM file..." >&2
        samtools index -@ ${THREADS} "${BAM_FILE}"
    fi
done

# --- 5. Variant Calling (lofreq) ---

for sample in "${SAMPLES[@]}"; do
    BAM_FILE="${RESULTS_DIR}/${sample}.bam"
    VCF_OUT="${RESULTS_DIR}/${sample}.vcf"
    
    if [[ -f "${BAM_FILE}" ]]; then
        echo "Calling variants for ${sample}..." >&2
        # lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o results/{sample}.vcf results/{sample}.bam
        lofreq call-parallel --pp-threads ${THREADS} -f "${REF_FA}" -o "${RESULTS_DIR}/${sample}.vcf" "${BAM_FILE}"
    else
        echo "Skipping variant calling for ${sample}: BAM file not found." >&2
    fi
done

# --- 6. VCF Compression and Indexing ---

for sample in "${SAMPLES[@]}"; do
    VCF_UNCOMPRESSED="${RESULTS_DIR}/${sample}.vcf"
    VCF_GZ="${RESULTS_DIR}/${sample}.vcf.gz"
    TBI_FILE="${RESULTS_DIR}/${sample}.vcf.gz.tbi"

    if [[ -f "${VCF_UNCOMPRESSED}" ]]; then
        echo "Compressing and indexing ${sample} VCF..." >&2
        bgzip -c "${VCF_UNCOMPRESSED}" > "${VCF_GZ}"
        tabix -p vcf "${VCF_GZ}"
        rm "${VCF_UNCOMPRESSED}"
    fi
done

# --- 7. Collapse Step ---

echo "Generating collapsed summary table..." >&2
HEADER="sample\tchrom\tpos\tref\talt\taf"
TEMP_OUTPUT=""

# Build the content for the collapsed TSV
for sample in "${SAMPLES[@]}"; do
    VCF_GZ="${RESULTS_DIR}/${sample}.vcf.gz"
    
    if [[ -f "${VCF_GZ}" ]]; then
        # bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/{sample}.vcf.gz
        # The sample name is prepended as a literal string in the format specifier.
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "${VCF_GZ}"
    fi
done | (
    # Prepend header to the combined stream
    echo -e "${HEADER}"
) > "${COLLAPSED_FILE}"

# Final check to ensure the collapsed file was written
if [[ -s "${COLLAPSED_FILE}" ]]; then
    echo "Workflow complete. Results available in ${RESULTS_DIR}/" >&2
else
    echo "Warning: Collapsed file ${COLLAPSED_FILE} is empty or missing." >&2
fi

exit 0