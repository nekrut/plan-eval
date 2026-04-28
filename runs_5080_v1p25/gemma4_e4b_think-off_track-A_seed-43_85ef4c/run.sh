#!/usr/bin/env bash
set -euo pipefail

# --- Configuration ---
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF_FA="data/ref/chrM.fa"
RAW_DIR="data/raw"
RES_DIR="results"
COLLAPSED_FILE="${RES_DIR}/collapsed.tsv"

# --- Setup and Idempotency Check ---

# Create results directory if it doesn't exist
mkdir -p "${RES_DIR}"

# Check if all required outputs exist to determine if work needs to be done.
# We check for the most complex/final artifact: the index of the collapsed file.
if [[ -f "${RES_DIR}/collapsed.tsv" ]] && \
   [[ -f "${RES_DIR}/M117-bl.vcf.gz.tbi" ]] && \
   [[ -f "${RES_DIR}/M117-ch.vcf.gz.tbi" ]] && \
   [[ -f "${RES_DIR}/M117C1-bl.vcf.gz.tbi" ]] && \
   [[ -f "${RES_DIR}/M117C1-ch.vcf.gz.tbi" ]]; then
    # If all primary artifacts exist, assume success and exit 0.
    exit 0
fi

# --- 1. Reference Indexing (Idempotent) ---

# 1a. samtools faidx
FAIDX_OUT="${RES_DIR}/chrM.fa.fai"
if [[ ! -f "${FAIDX_OUT}" ]] || find "${RES_DIR}" -maxdepth 1 -name "*.fai" -print | grep -q "chrM.fa.fai"; then
    samtools faidx "${REF_FA}"
    echo "Indexed reference FASTA file."
fi

# 1b. bwa index
BWA_INDEX_FILES=("${RES_DIR}/chrM.fa.amb" "${RES_DIR}/chrM.fa.ann" "${RES_DIR}/chrM.fa.bwt" "${RES_DIR}/chrM.fa.pac" "${RES_DIR}/chrM.fa.sa")
if ! ls "${RES_DIR}"/chrM.fa.* 1> /dev/null 2>&1; then
    bwa index "${REF_FA}"
    echo "Indexed reference using bwa."
fi

# --- 2. Per-sample Alignment, Sorting, and Indexing ---
for sample in "${SAMPLES[@]}"; do
    echo "Processing sample: ${sample}"

    # Define input files
    R1="${RAW_DIR}/${sample}_1.fq.gz"
    R2="${RAW_DIR}/${sample}_2.fq.gz"
    BAM_OUT="${RES_DIR}/${sample}.bam"
    BAI_OUT="${RES_DIR}/${sample}.bam.bai"

    # Check if BAM exists and is recent enough (simple check for existence)
    if [[ -f "${BAM_OUT}" ]] && [[ -f "${BAI_OUT}" ]]; then
        # A more robust check would compare modification times, but for simplicity and idempotency,
        # we rely on the fact that if the VCF/TBI exists, the BAM/BAI should too.
        : # Skip if BAM/BAI exist
    else
        # 2a. Alignment (bwa mem)
        # Read Group format: -R "@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA"
        RG_LINE="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
        
        # Use printf to ensure the literal backslash-t sequence is passed correctly to bwa
        bwa mem -t ${THREADS} -R "${RG_LINE}" "${REF_FA}" "${R1}" "${R2}" | \
        # 2b. SAM -> Sorted BAM
        samtools sort -@ ${THREADS} -o "${BAM_OUT}" -
        
        # 2c. BAM Indexing
        samtools index -@ ${THREADS} "${BAM_OUT}"
        echo "Finished alignment and indexing for ${sample}."
    fi
done

# --- 3. Variant Calling (lofreq) ---
for sample in "${SAMPLES[@]}"; do
    VCF_OUT="${RES_DIR}/${sample}.vcf"
    VCFGZ_OUT="${RES_DIR}/${sample}.vcf.gz"
    TBI_OUT="${RES_DIR}/${sample}.vcf.gz.tbi"
    
    # Check if VCF.gz.tbi exists, implying the whole pipeline for this sample ran.
    if [[ -f "${TBI_OUT}" ]]; then
        : # Skip
    else
        # 3a. Variant Calling
        lofreq call-parallel --pp-threads ${THREADS} -f "${REF_FA}" -o "${RES_DIR}/${sample}.vcf" "${RES_DIR}/${sample}.bam"
        
        # 3b. Compression and Indexing
        bgzip -c "${RES_DIR}/${sample}.vcf" > "${VCFGZ_OUT}"
        tabix -p vcf "${VCFGZ_OUT}"
        
        # 3c. Cleanup intermediate file
        rm "${RES_DIR}/${sample}.vcf"
        echo "Finished variant calling and indexing for ${sample}."
    fi
done

# --- 4. Collapse Step (results/collapsed.tsv) ---

# Check if the collapsed file exists and if any input VCF is newer than it.
# We check the modification time of the newest input VCF against the output.
INPUT_VCFS=("$RES_DIR/${SAMPLES[0]}.vcf.gz" "${RES_DIR}/${SAMPLES[1]}.vcf.gz" "${RES_DIR}/${SAMPLES[2]}.vcf.gz" "${RES_DIR}/${SAMPLES[3]}.vcf.gz")
NEWEST_INPUT_TIME=0
for vcf in "${INPUT_VCFS[@]}"; do
    if [[ -f "$vcf" ]]; then
        # Get modification time in seconds since epoch
        mtime=$(stat -c %Y "$vcf")
        if (( mtime > NEWEST_INPUT_TIME )); then
            NEWEST_INPUT_TIME=$mtime
        fi
    fi
done

if [[ -f "${COLLAPSED_FILE}" ]] && [[ $(stat -c %Y "${COLLAPSED_FILE}") -ge "${NEWEST_INPUT_TIME:-0}" ]]; then
    : # Skip, output is up to date
else
    echo "Generating collapsed variant table: ${COLLAPSED_FILE}"
    
    # Header
    {
        echo -e "sample\tchrom\tpos\tref\talt\taf"
        
        # Process each sample
        for sample in "${SAMPLES[@]}"; do
            VCF_FILE="${RES_DIR}/${sample}.vcf.gz"
            
            if [[ -f "${VCF_FILE}" ]]; then
                # bcftools query format: sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF
                bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "${VCF_FILE}"
            fi
        done
    } | sort -k 2,2 -k 1,1 > "${COLLAPSED_FILE}"
    
    echo "Successfully created ${COLLAPSED_FILE}"
fi

exit 0