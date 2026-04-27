#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF_FA="data/ref/chrM.fa"
RESULTS_DIR="results"
COLLAPSED_FILE="${RESULTS_DIR}/collapsed.tsv"

# --- 1. Setup and Idempotency Check ---

mkdir -p "${RESULTS_DIR}"

# Check if all required outputs exist. If they do, we assume success and exit 0.
# This is a simplified check; a full check would compare timestamps, but for simplicity
# and robustness against partial runs, checking for the final collapsed file is sufficient
# to satisfy the "exit 0 if populated" requirement, provided the steps are ordered.
if [[ -f "${COLLAPSED_FILE}" ]] && [[ -f "${REF_FA}.fai" ]] && [[ -f "${REF_FA}.bwt" ]]; then
    # Check if the most recent input files are older than the output, or if the output is missing.
    # If the final file exists, we assume the process was successful enough to exit cleanly.
    # A more rigorous check would compare all timestamps, but this satisfies the requirement.
    echo "Results directory appears populated. Exiting successfully." >&2
    exit 0
fi

# --- 2. Reference Indexing ---

echo "Indexing reference genome..." >&2
if [[ ! -f "${REF_FA}.fai" ]]; then
    samtools faidx "${REF_FA}"
    if [ $? -ne 0 ]; then
        echo "Error indexing reference FASTA file. Aborting." >&2
        exit 1
    fi
fi

# Check for BWA index files (bwt, sa, etc.)
if [[ ! -f "${REF_FA}.bwt" ]]; then
    bwa index "${REF_FA}"
    if [ $? -ne 0 ]; then
        echo "Error indexing reference with BWA. Aborting." >&2
        exit 1
    fi
fi

# --- 3. & 4. Per-sample Alignment and Sorting ---

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: ${SAMPLE}" >&2
    
    R1="${data/raw/${SAMPLE}_1.fq.gz"
    R2="${data/raw/${SAMPLE}_2.fq.gz"
    BAM_OUT="${RESULTS_DIR}/${SAMPLE}.bam"
    
    # Check if BAM exists and is up-to-date (simple check: if it exists, assume done)
    if [[ -f "${BAM_OUT}" ]]; then
        echo "  Skipping alignment for ${SAMPLE}: ${BAM_OUT} already exists." >&2
        continue
    fi

    # Construct the Read Group header string carefully to avoid shell expansion issues
    # The required format is: -R "@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA"
    # We must pass the literal backslash-t sequence.
    RG_HEADER="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA"

    echo "  Aligning reads for ${SAMPLE}..." >&2
    
    # bwa mem output is piped directly to samtools sort
    bwa mem -t ${THREADS} -R "${RG_HEADER}" "${REF_FA}" "${R1}" "${R2}" | \
    samtools sort -@ ${THREADS} -o "${BAM_OUT}" -T "${RESULTS_DIR}/${SAMPLE}.tmp.sorted.bam" >/dev/null
    
    # Rename the temporary sorted file to the final expected name
    mv "${RESULTS_DIR}/${SAMPLE}.tmp.sorted.bam" "${BAM_OUT}"
    
    if [ $? -ne 0 ]; then
        echo "Error during alignment/sorting for ${SAMPLE}. Aborting." >&2
        exit 1
    fi
done

# --- 5. BAM Indexing ---

for SAMPLE in "${SAMPLES[@]}"; do
    BAM_OUT="${RESULTS_DIR}/${SAMPLE}.bam"
    BAI_OUT="${RESULTS_DIR}/${SAMPLE}.bam.bai"
    
    if [[ -f "${BAI_OUT}" ]]; then
        echo "  Skipping indexing for ${SAMPLE}: ${BAI_OUT} already exists." >&2
        continue
    fi
    
    echo "  Indexing BAM for ${SAMPLE}..." >&2
    samtools index -@ ${THREADS} "${BAM_OUT}"
done

# --- 6. Variant Calling with lofreq ---

for SAMPLE in "${SAMPLES[@]}"; do
    BAM_OUT="${RESULTS_DIR}/${SAMPLE}.bam"
    VCF_OUT_GZ="${RESULTS_DIR}/${SAMPLE}.vcf.gz"
    VCF_TEMP="${RESULTS_DIR}/${SAMPLE}.vcf"
    TBI_OUT="${RESULTS_DIR}/${SAMPLE}.vcf.gz.tbi"
    
    # Check if VCF.gz exists and is indexed
    if [[ -f "${TBI_OUT}" ]]; then
        echo "  Skipping variant calling for ${SAMPLE}: ${VCF_OUT_GZ} and ${TBI_OUT} exist." >&2
        continue
    fi

    echo "  Calling variants for ${SAMPLE}..." >&2
    
    # 6. Call variants (outputs uncompressed VCF to temp file)
    lofreq call-parallel --threads ${THREADS} \
        --reference "${REF_FA}" \
        --bam "${BAM_OUT}" \
        --output "${VCF_TEMP}"
    
    if [ $? -ne 0 ]; then
        echo "Error during lofreq calling for ${SAMPLE}. Aborting." >&2
        exit 1
    fi
    
    # 7. VCF compression and indexing
    echo "  Compressing and indexing VCF for ${SAMPLE}..." >&2
    bgzip -c "${VCF_TEMP}" > "${VCF_OUT_GZ}"
    tabix -p vcf "${VCF_OUT_GZ}"
    
    # Cleanup intermediate file
    rm "${VCF_TEMP}"
done

# --- 8. Collapse Step ---

# Check if the final collapsed file exists and is up-to-date
if [[ -f "${COLLAPSED_FILE}" ]]; then
    # Check if any input VCF.gz is newer than the collapsed file
    local needs_rebuild=0
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF_INPUT="${RESULTS_DIR}/${SAMPLE}.vcf.gz"
        if [[ ! -f "${VCF_INPUT}" ]] || [[ $(stat -c %Y "${VCF_INPUT}") -gt $(stat -c %Y "${COLLAPSED_FILE}") ]]; then
            needs_rebuild=1
            break
        fi
    done
    if [[ $needs_rebuild -eq 0 ]]; then
        echo "Collapsed table ${COLLAPSED_FILE} is up-to-date. Skipping." >&2
        exit 0
    fi
fi

echo "  Generating collapsed variant table..." >&2

# Header
HEADER="sample\tchrom\tpos\tref\talt\taf"

# Build the query command structure
QUERY_CMDS=""
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_INPUT="${RESULTS_DIR}/${SAMPLE}.vcf.gz"
    # Format string: {sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n
    # The sample name is prepended literally to the output stream.
    QUERY_CMDS+=( "bcftools query -f \"${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n\" ${VCF_INPUT}" )
done

# Execute all queries and concatenate, prepending header
{
    echo -e "${HEADER}"
    for cmd in "${QUERY_CMDS[@]}"; do
        eval "$cmd"
    done
} | sort -k1,1 -k2,2 -k3,3 > "${COLLAPSED_FILE}"

if [ $? -eq 0 ]; then
    echo "Successfully created ${COLLAPSED_FILE}" >&2
else
    echo "Error creating ${COLLAPSED_FILE}" >&2
    exit 1
fi

exit 0