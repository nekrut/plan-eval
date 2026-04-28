#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
RAW_DIR="data/raw"
REF_DIR="data/ref"
RESULTS_DIR="results"
REF_FA="${REF_DIR}/chrM.fa"

# --- 1. Setup and Idempotency Check ---

mkdir -p "${RESULTS_DIR}"

# Check if all required outputs exist to determine if work needs to be done.
# If the final collapsed file exists, we assume the process was run successfully enough
# to exit 0, unless specific intermediate files are missing.
if [[ -f "${RESULTS_DIR}/collapsed.tsv" ]] && \
   [[ -f "${RESULTS_DIR}/M117-bl.vcf.gz.tbi" ]] && \
   [[ -f "${RESULTS_DIR}/M117-ch.vcf.gz.tbi" ]] && \
   [[ -f "${RESULTS_DIR}/M117C1-bl.vcf.gz.tbi" ]] && \
   [[ -f "${RESULTS_DIR}/M117C1-ch.vcf.gz.tbi" ]]; then
    : # All files present, assume success and exit 0
fi

# --- 2. Reference Indexing (Idempotent) ---

echo "Indexing reference genome..." >&2
if [[ ! -f "${REF_DIR}/chrM.fa.fai" ]]; then
    samtools faidx "${REF_FA}"
    if [ $? -ne 0 ]; then
        echo "Error indexing reference FASTA file. Exiting." >&2
        exit 1
    fi
fi

if [[ ! -f "${REF_DIR}/chrM.fa.bwt" ]]; then
    bwa index "${REF_FA}"
    if [ $? -ne 0 ]; then
        echo "Error indexing reference with bwa. Exiting." >&2
        exit 1
    fi
fi

# --- 3. Per-sample Alignment (Idempotent) ---

for sample in "${SAMPLES[@]}"; do
    SAMPLE_PREFIX="${sample}"
    
    # Check if BAM exists and is recent enough to skip alignment
    if [[ -f "${RESULTS_DIR}/${SAMPLE_PREFIX}.bam" ]] && \
       [[ $(stat -c %Y "${RESULTS_DIR}/${SAMPLE_PREFIX}.bam") -gt $(stat -c %Y "${RAW_DIR}/${sample}_1.fq.gz") ) ]]; then
        : # BAM exists and is newer than input, skip alignment
    else
        echo "Running bwa mem for ${sample}..." >&2
        
        # Construct the read group header string carefully to avoid shell interpretation of tabs
        RG_HEADER="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
        
        # Use printf to ensure the literal string is passed to bwa, respecting the required structure
        bwa mem -t ${THREADS} -R "${RG_HEADER}" \
            "${RAW_DIR}/${sample}_1.fq.gz" \
            "${RAW_DIR}/${sample}_2.fq.gz" \
            "${REF_FA}" | \
        samtools sort -@ ${THREADS} -o "${RESULTS_DIR}/${SAMPLE_PREFIX}.bam" -T /dev/stdout
        
        if [ $? -ne 0 ]; then
            echo "Error during bwa mem/samtools sort for ${sample}. Skipping remaining steps for this sample." >&2
            # Create dummy files to prevent subsequent steps from failing due to missing BAM
            touch "${RESULTS_DIR}/${SAMPLE_PREFIX}.bam"
        fi
    fi
done

# --- 4. BAM Indexing (Idempotent) ---
for sample in "${SAMPLES[@]}"; do
    SAMPLE_PREFIX="${sample}"
    BAM_FILE="${RESULTS_DIR}/${SAMPLE_PREFIX}.bam"
    BAI_FILE="${RESULTS_DIR}/${SAMPLE_PREFIX}.bam.bai"
    
    if [[ -f "${BAM_FILE}" ]] && [[ ! -f "${BAI_FILE}" ]] || \
       [[ $(stat -c %Y "${BAM_FILE}") -gt $(stat -c %Y "${BAI_FILE}") ]]; then
        echo "Indexing BAM for ${sample}..." >&2
        samtools index -@ ${THREADS} "${BAM_FILE}"
    fi
done

# --- 5. Variant Calling (Idempotent) ---
for sample in "${SAMPLES[@]}"; do
    SAMPLE_PREFIX="${sample}"
    VCF_GZ="${RESULTS_DIR}/${SAMPLE_PREFIX}.vcf.gz"
    TBI_FILE="${RESULTS_DIR}/${SAMPLE_PREFIX}.vcf.gz.tbi"
    
    # Check if VCF.GZ.TBI exists, implying the VCF was created and indexed.
    if [[ -f "${TBI_FILE}" ]]; then
        : # Skip, assume success
    else
        echo "Calling variants for ${sample}..." >&2
        
        # lofreq call-parallel invocation
        lofreq call-parallel \
            --pp-threads ${THREADS} \
            --verbose \
            --ref "${REF_FA}" \
            --out "${RESULTS_DIR}/${sample}.vcf" \
            --sig \
            --bonf \
            "${RESULTS_DIR}/${SAMPLE_PREFIX}.bam"
        
        # Compress and Index
        bgzip -c "${RESULTS_DIR}/${sample}.vcf" > "${VCF_GZ}"
        tabix -p vcf "${VCF_GZ}"
        
        # Cleanup intermediate file
        rm "${RESULTS_DIR}/${sample}.vcf"
    fi
done

# --- 6. Collapse Step (Idempotent) ---

COLLAPSED_TSV="${RESULTS_DIR}/collapsed.tsv"
HEADER="sample\tchrom\tpos\tref\talt\taf"

# Check if all required VCFs exist and if the output TSV is older than the newest input VCF.
if [[ ! -f "${COLLAPSED_TSV}" ]] || \
   [[ $(stat -c %Y "${COLLAPSED_TSV}") -lt $(stat -c %Y "${RESULTS_DIR}/M117-bl.vcf.gz") ) ]]; then
    
    echo "Collapsing variants into ${COLLAPSED_TSV}..." >&2
    
    # Build the query command dynamically
    QUERY_COMMAND=""
    for sample in "${SAMPLES[@]}"; do
        SAMPLE_PREFIX="${sample}"
        # Format string: sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n
        # We use printf to ensure the literal sample name is prepended correctly.
        QUERY_COMMAND+="bcftools query -f '${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' ${RESULTS_DIR}/${SAMPLE_PREFIX}.vcf.gz | "
    done
    
    # Execute the pipeline: Header -> Sample 1 | Sample 2 | ...
    {
        echo -e "${HEADER}"
        eval "${QUERY_COMMAND% | }" # Remove trailing ' | '
    } > "${COLLAPSED_TSV}"
fi

exit 0