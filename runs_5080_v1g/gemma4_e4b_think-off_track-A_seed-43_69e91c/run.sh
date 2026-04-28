#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF_FA="data/ref/chrM.fa"
RAW_DIR="data/raw"
RESULTS_DIR="results"
COLLAPSED_FILE="${RESULTS_DIR}/collapsed.tsv"

# --- 1. Setup and Idempotency Check ---

mkdir -p "${RESULTS_DIR}"

# Check if all required outputs exist to determine if work needs to be done.
# If the final collapsed file exists, we assume the process was run successfully enough
# to exit 0, unless a critical input is missing.

# --- 2. Reference Indexing ---

if [ ! -f "${REF_FA}.fai" ] || [ ! -f "${REF_FA}.bwt" ]; then
    # samtools faidx
    if ! samtools faidx "${REF_FA}"; then
        echo "Error indexing reference FASTA file with samtools faidx. Exiting." >&2
        exit 1
    fi
    # bwa index
    if ! bwa index "${REF_FA}"; then
        echo "Error indexing reference FASTA file with bwa index. Exiting." >&2
        exit 1
    fi
fi

# --- 3. Per-sample Alignment (bwa mem) ---

for sample in "${SAMPLES[@]}"; do
    SAMPLE_PREFIX="${sample}"
    R1="${RAW_DIR}/${SAMPLE_PREFIX}_1.fq.gz"
    R2="${RAW_DIR}/${SAMPLE_PREFIX}_2.fq.gz"
    BAM_OUT="${RESULTS_DIR}/${SAMPLE_PREFIX}.bam"

    if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
        echo "Skipping ${sample}: Input FASTQ files not found." >&2
        continue
    fi

    # Check if BAM exists and is recent enough (simple check for existence)
    if [ -f "${BAM_OUT}" ] && [ $(stat -c %Y "${BAM_OUT}") -gt 0 ]; then
        # Assume success if BAM exists, skip alignment
        :
    else
        # Read Group construction: -R "@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA"
        # Note: The literal backslash-t must be passed correctly to bwa.
        RG_LINE="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
        
        # Pipe bwa mem output directly to samtools sort
        if ! bwa mem -t ${THREADS} -R "${RG_LINE}" "${REF_FA}" "${R1}" "${R2}" | \
            samtools sort -@ ${THREADS} -o "${BAM_OUT}"; then
            echo "Error during bwa mem and sorting for ${sample}." >&2
            exit 1
        fi
    fi
done

# --- 4. BAM Indexing ---

for sample in "${SAMPLES[@]}"; do
    SAMPLE_PREFIX="${sample}"
    BAM_OUT="${RESULTS_DIR}/${SAMPLE_PREFIX}.bam"
    BAI_OUT="${RESULTS_DIR}/${SAMPLE_PREFIX}.bam.bai"

    if [ -f "${BAM_OUT}" ] && [ ! -f "${BAI_OUT}" ]; then
        if ! samtools index -@ ${THREADS} "${BAM_OUT}"; then
            echo "Error indexing BAM file for ${sample}." >&2
            exit 1
        fi
    fi
done

# --- 5. Variant Calling (lofreq) ---

for sample in "${SAMPLES[@]}"; do
    SAMPLE_PREFIX="${sample}"
    BAM_IN="${RESULTS_DIR}/${SAMPLE_PREFIX}.bam"
    VCF_OUT="${RESULTS_DIR}/${SAMPLE_PREFIX}.vcf.gz"
    TBI_OUT="${RESULTS_DIR}/${SAMPLE_PREFIX}.vcf.gz.tbi"

    if [ ! -f "${BAM_IN}" ]; then
        echo "Skipping variant calling for ${sample}: BAM file not found." >&2
        continue
    fi

    # Check if VCF.gz and TBI exist
    if [ -f "${VCF_OUT}" ] && [ -f "${TBI_OUT}" ]; then
        : # Skip
    else
        # lofreq call-parallel invocation
        if ! lofreq call-parallel \
            --pp-threads ${THREADS} \
            --verbose \
            --ref "${REF_FA}" \
            --out "${RESULTS_DIR}/${sample}.vcf" \
            --sig \
            --bonf \
            "${BAM_IN}"; then
            echo "Error during lofreq call-parallel for ${sample}." >&2
            exit 1
        fi
        
        # Compress and index
        if ! bgzip -c "${RESULTS_DIR}/${sample}.vcf" > "${VCF_OUT}"; then
            echo "Error compressing VCF for ${sample}." >&2
            exit 1
        fi
        
        if ! tabix -p vcf "${VCF_OUT}"; then
            echo "Error indexing VCF for ${sample}." >&2
            exit 1
        fi
        
        # Cleanup intermediate file
        rm "${RESULTS_DIR}/${sample}.vcf"
    fi
done

# --- 6. Collapse Step ---

# Check if all VCFs exist before attempting collapse
ALL_VCFS_PRESENT=true
for sample in "${SAMPLES[@]}"; do
    if [ ! -f "${RESULTS_DIR}/${sample}.vcf.gz" ]; then
        ALL_VCFS_PRESENT=false
        break
    fi
done

if ! $ALL_VCFS_PRESENT; then
    echo "Skipping collapse step: One or more sample VCFs are missing." >&2
    exit 0
fi

# Check if the final output file is outdated
if [ -f "${COLLAPSED_FILE}" ] && [ $(stat -c %Y "${COLLAPSED_FILE}") -gt 0 ]; then
    # Check if any input VCF is newer than the output
    NEWER_INPUT=false
    for sample in "${SAMPLES[@]}"; do
        if [ -f "${RESULTS_DIR}/${sample}.vcf.gz" ] && [ $(stat -c %Y "${RESULTS_DIR}/${sample}.vcf.gz}") -gt $(stat -c %Y "${COLLAPSED_FILE}") ]; then
            NEWER_INPUT=true
            break
        fi
    done
    if ! $NEWER_INPUT; then
        : # Skip collapse
    else
        : # Proceed with collapse
    fi
else
    : # Proceed with collapse (file missing or zero size)
fi


# Execute collapse if necessary
{
    # Header
    echo -e "sample\tchrom\tpos\tref\talt\taf"
    
    # Data lines
    for sample in "${SAMPLES[@]}"; do
        VCF_IN="${RESULTS_DIR}/${sample}.vcf.gz"
        # bcftools query format: {sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "${VCF_IN}"
    done
} | sort -k1,1 -k2,2 -k3,3 > "${COLLAPSED_FILE}"

# Final exit 0 ensures idempotency on success
exit 0