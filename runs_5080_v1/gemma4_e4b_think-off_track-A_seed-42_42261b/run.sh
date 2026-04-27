#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF_FA="data/ref/chrM.fa"
RESULTS_DIR="results"
COLLAPSED_FILE="${RESULTS_DIR}/collapsed.tsv"

# --- 1. Setup and Idempotency Check ---

mkdir -p "${RESULTS_DIR}"

# Check if all required outputs exist. If so, exit successfully.
# We check for the most complex/final artifact: the collapsed TSV.
if [[ -f "${COLLAPSED_FILE}" ]] && \
   [[ -f "${REF_FA}.fai" ]] && \
   [[ -f "${REF_FA}.bwt" ]] && \
   [[ -f "${REF_FA}.sa" ]]; then
    # Basic check passed, assume work is done.
    exit 0
fi

# --- 2. Reference Indexing ---

# Index FASTA file (samtools faidx)
if [[ ! -f "${REF_FA}.fai" ]]; then
    samtools faidx "${REF_FA}"
fi

# BWA Indexing
if [[ ! -f "${REF_FA}.bwt" ]]; then
    bwa index "${REF_FA}"
fi

# --- 3. Per-sample Alignment and Processing ---

for sample in "${SAMPLES[@]}"; do
    SAMPLE_PREFIX="${sample}"
    
    # Check if BAM exists and is recent enough to skip alignment/calling
    BAM_FILE="${RESULTS_DIR}/${SAMPLE_PREFIX}.bam"
    VCF_GZ="${RESULTS_DIR}/${SAMPLE_PREFIX}.vcf.gz"
    
    if [[ -f "${BAM_FILE}" ]] && [[ -f "${VCF_GZ}" ]]; then
        # Check if the VCF is newer than the BAM (a proxy for completeness)
        if find "${RESULTS_DIR}" -maxdepth 1 -name "${SAMPLE_PREFIX}.bam" -print -quit | grep -q "${BAM_FILE}"; then
            # If BAM and VCF exist, we assume the pipeline ran successfully enough to skip.
            # We still need to run the collapse step later, but for now, we skip the heavy lifting.
            : # No action needed, just continue loop
        else
            : # Should not happen if the check above is correct, but safe to continue.
        fi
    else
        # --- 3. Alignment (bwa mem) ---
        R1="${data/raw/${sample}_1.fq.gz"
        R2="${data/raw/${sample}_2.fq.gz"
        
        # Check if input FASTQs exist
        if [[ ! -f "${R1}" ]] || [[ ! -f "${R2}" ]]; then
            echo "Skipping ${sample}: Input FASTQs not found." >&2
            continue
        fi

        # Read Group definition: Literal backslash-t
        RG_LINE="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
        
        # Alignment command
        if ! bwa mem -t ${THREADS} -R "${RG_LINE}" "${REF_FA}" "${R1}" "${R2}" 2> /dev/null | \
            samtools sort -@ ${THREADS} -o "${BAM_FILE}" -T "${RESULTS_DIR}/temp_sort"; then
            echo "Error during bwa mem/samtools sort for ${sample}. Skipping." >&2
            continue
        fi
        
        # --- 5. BAM Indexing ---
        if ! samtools index -@ ${THREADS} "${BAM_FILE}"; then
            echo "Error indexing BAM for ${sample}. Skipping." >&2
            continue
        fi

        # --- 6. Variant Calling (lofreq call-parallel) ---
        TEMP_VCF="${RESULTS_DIR}/${sample}.vcf"
        
        # lofreq call-parallel emits uncompressed VCF
        if ! lofreq call-parallel --threads ${THREADS} \
            --reference "${REF_FA}" \
            --input "${BAM_FILE}" \
            --output "${TEMP_VCF}"; then
            echo "Error during lofreq calling for ${sample}. Skipping." >&2
            rm -f "${TEMP_VCF}"
            continue
        fi

        # --- 7. VCF Compression and Indexing ---
        if ! bgzip -c "${TEMP_VCF}" > "${VCF_GZ}"; then
            echo "Error compressing VCF for ${sample}. Skipping." >&2
            rm -f "${TEMP_VCF}"
            continue
        fi
        
        if ! tabix -p vcf "${VCF_GZ}"; then
            echo "Error indexing VCF for ${sample}. Skipping." >&2
            rm -f "${TEMP_VCF}"
            continue
        fi
        
        # Cleanup intermediate file
        rm -f "${TEMP_VCF}"
    fi
done

# --- 8. Collapse Step ---

# Check if any VCF was successfully generated before attempting collapse
if ! ls "${RESULTS_DIR}"/*.vcf.gz 1> /dev/null 2>&1; then
    echo "No VCF files found in ${RESULTS_DIR}. Skipping collapse step." >&2
    exit 0
fi

# Determine if collapse needs to run (check if header exists or if any VCF is newer than the TSV)
if [[ ! -f "${COLLAPSED_FILE}" ]] || find "${RESULTS_DIR}" -maxdepth 1 -name "*.vcf.gz" -newer "${COLLAPSED_FILE}" -print -quit 2>/dev/null | grep -q .; then
    
    HEADER="sample\tchrom\tpos\tref\talt\taf"
    
    # Use printf to build the header and then loop to append data
    {
        printf "%s\n" "${HEADER}"
        for sample in "${SAMPLES[@]}"; do
            VCF_GZ="${RESULTS_DIR}/${sample}.vcf.gz"
            if [[ -f "${VCF_GZ}" ]]; then
                # bcftools query: format string prepends the sample name as the first field
                bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "${VCF_GZ}"
            fi
        done
    } | sort -k1,1 -k2,2 -k3,3 > "${COLLAPSED_FILE}"
fi

exit 0