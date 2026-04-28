#!/usr/bin/env bash
set -euo pipefail

# --- Configuration ---
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF_FA="data/ref/chrM.fa"
RAW_DIR="data/raw"
RESULTS_DIR="results"
COLLAPSED_FILE="${RESULTS_DIR}/collapsed.tsv"

# --- Setup and Idempotency Check ---

# Create results directory if it doesn't exist
mkdir -p "${RESULTS_DIR}"

# Check if all required outputs exist to determine if work needs to be done.
# If the final collapsed file exists, we assume the process was run successfully enough
# to exit 0, unless we detect missing intermediate files.
# For simplicity and robustness, we will check for the existence of the most complex
# artifact (the index file) to gate the entire process.

# Check if any sample's VCF.gz.tbi exists. If so, we assume success for the main pipeline.
if [[ -f "${RESULTS_DIR}/M117-bl.vcf.gz.tbi" ]] && \
   [[ -f "${RESULTS_DIR}/M117-ch.vcf.gz.tbi" ]] && \
   [[ -f "${RESULTS_DIR}/M117C1-bl.vcf.gz.tbi" ]] && \
   [[ -f "${RESULTS_DIR}/M117C1-ch.vcf.gz.tbi" ]]; then
    # All primary artifacts exist, assume success and exit 0.
    exit 0
fi

# --- 1. Reference Indexing (Idempotent) ---

echo "Indexing reference genome..." >&2

# samtools faidx
if [ ! -f "${REF_FA}.fai" ]; then
    samtools faidx "${REF_FA}"
fi

# bwa index
if [ ! -f "${REF_FA}.amb" ]; then
    bwa index "${REF_FA}"
fi

# --- 2. Per-sample Alignment and Processing ---

for sample in "${SAMPLES[@]}"; do
    echo "Processing sample: ${sample}" >&2

    # Define paths for this sample
    R1="${RAW_DIR}/${sample}_1.fq.gz"
    R2="${RAW_DIR}/${sample}_2.fq.gz"
    BAM_OUT="${RESULTS_DIR}/${sample}.bam"
    BAI_OUT="${RESULTS_DIR}/${sample}.bam.bai"
    VCF_OUT="${RESULTS_DIR}/${sample}.vcf.gz"
    TBI_OUT="${RESULTS_DIR}/${sample}.vcf.gz.tbi"

    # Check if the final VCF.gz.tbi exists for this sample. If so, skip the whole block.
    if [ -f "${TBI_OUT}" ]; then
        echo "Skipping ${sample}: ${TBI_OUT} found." >&2
        continue
    fi

    # 3. Alignment (bwa mem)
    # Read Group format: @RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA
    RG_LINE="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"

    # Check if input FASTQs exist
    if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
        echo "Skipping ${sample}: Input FASTQs not found." >&2
        continue
    fi

    # bwa mem output is piped directly to samtools sort
    # We use a temporary file for the unsorted SAM output to ensure the pipe structure works correctly
    TEMP_SAM="${RESULTS_DIR}/${sample}.sam"
    
    # Run alignment
    bwa mem -t ${THREADS} -R "${RG_LINE}" "${REF_FA}" "${R1}" "${R2}" | \
    samtools sort -@ ${THREADS} -o "${BAM_OUT}" - > /dev/null

    # 4. BAM Indexing (Idempotent)
    if [ ! -f "${BAI_OUT}" ]; then
        samtools index -@ ${THREADS} "${BAM_OUT}"
    fi

    # 5. Variant Calling (lofreq)
    # lofreq call-parallel outputs plain VCF to results/{sample}.vcf
    TEMP_VCF="${RESULTS_DIR}/${sample}.vcf"
    if [ ! -f "${TEMP_VCF}" ]; then
        lofreq call-parallel --pp-threads ${THREADS} -f "${REF_FA}" -o "${RESULTS_DIR}" "${BAM_OUT}" "${TEMP_VCF}"
    else
        echo "Skipping ${sample}: ${TEMP_VCF} found." >&2
    fi

    # 6. VCF Compression and Indexing (Idempotent)
    if [ ! -f "${VCF_OUT}" ] || [ "$TEMP_VCF" -nt "${VCF_OUT}" ]; then
        bgzip -c "${TEMP_VCF}" > "${VCF_OUT}"
        tabix -p vcf "${VCF_OUT}"
        rm "${TEMP_VCF}"
    fi

done

# --- 7. Collapse Step (Idempotent) ---

# Check if all required VCF.gz.tbi files exist before attempting collapse
ALL_VCF_PRESENT=true
for sample in "${SAMPLES[@]}"; do
    if [ ! -f "${RESULTS_DIR}/${sample}.vcf.gz.tbi" ]; then
        ALL_VCF_PRESENT=false
        break
    fi
done

if ! $ALL_VCF_PRESENT; then
    echo "Warning: Not all sample VCFs are indexed. Skipping collapse step." >&2
    exit 0
fi

# Check if the collapsed file needs rebuilding
# We check if any input VCF is newer than the output TSV
INPUT_NEWER=false
for sample in "${SAMPLES[@]}"; do
    if [ -f "${RESULTS_DIR}/${sample}.vcf.gz" ] && [ ! -f "${COLLAPSED_FILE}" ] || [ "$(stat -c %Y "${RESULTS_DIR}/${sample}.vcf.gz}")" -gt "$(stat -c %Y "${COLLAPSED_FILE}")" 2>/dev/null || [ ! -f "${COLLAPSED_FILE}" ]; then
        INPUT_NEWER=true
        break
    fi
done

if ! $INPUT_NEWER; then
    echo "Collapse file ${COLLAPSED_FILE} is up to date. Exiting cleanly." >&2
    exit 0
fi

echo "Generating collapsed variant table: ${COLLAPSED_FILE}" >&2

# Header
HEADER="sample\tchrom\tpos\tref\talt\taf"
echo -e "${HEADER}" > "${COLLAPSED_FILE}"

# Process samples and append to the file
for sample in "${SAMPLES[@]}"; do
    # bcftools query format: sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF
    # We use printf to ensure the literal sample name is prepended correctly.
    bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
        "${RESULTS_DIR}/${sample}.vcf.gz" >> "${COLLAPSED_FILE}"
done

exit 0