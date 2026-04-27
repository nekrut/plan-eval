#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF_FA="data/ref/chrM.fa"
RESULTS_DIR="results"
COLLAPSED_TSV="${RESULTS_DIR}/collapsed.tsv"

# --- 1. Setup and Idempotency Check ---

mkdir -p "${RESULTS_DIR}"

# Check if all required outputs exist. If they do, assume success and exit 0.
# This check is complex due to dependencies, so we rely on the explicit checks below
# to ensure no work is done if the final state is present.

# --- 2. Reference Indexing ---

echo "Indexing reference genome..." >&2
if [ ! -f "${REF_FA}.fai" ]; then
    samtools faidx "${REF_FA}"
    if [ $? -ne 0 ]; then
        echo "Error indexing reference FASTA file: ${REF_FA}" >&2
        exit 1
    fi
fi

if [ ! -f "${REF_FA}.bwt" ]; then
    bwa index "${REF_FA}"
    if [ $? -ne 0 ]; then
        echo "Error indexing reference with bwa." >&2
        exit 1
    fi
fi

# --- 3. Per-sample Alignment (bwa mem) ---

for sample in "${SAMPLES[@]}"; do
    echo "Processing sample: ${sample}" >&2
    R1="${data/raw/${sample}_1.fq.gz"
    R2="${data/raw/${sample}_2.fq.gz"
    
    # Check if inputs exist
    if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
        echo "Skipping ${sample}: Input FASTQ files not found." >&2
        continue
    fi

    # Check if final BAM exists and is recent enough (simple check for existence)
    BAM_OUT="${RESULTS_DIR}/${sample}.bam"
    if [ -f "${BAM_OUT}" ] && [ -s "${BAM_OUT}" ]; then
        echo "Skipping alignment for ${sample}: ${BAM_OUT} already exists." >&2
        continue
    fi

    # Read Group definition: literal backslash-t
    RG_LINE="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
    
    # Run bwa mem
    bwa mem -t ${THREADS} -R "${RG_LINE}" "${REF_FA}" "${R1}" "${R2}" | \
    samtools sort -@ ${THREADS} -o "${BAM_OUT}" -

    if [ $? -ne 0 ]; then
        echo "Error during bwa mem or samtools sort for ${sample}." >&2
        exit 1
    fi
done

# --- 4. BAM Indexing ---

for sample in "${SAMPLES[@]}"; do
    BAM_OUT="${RESULTS_DIR}/${sample}.bam"
    if [ -f "${BAM_OUT}" ] && [ -s "${BAM_OUT}" ]; then
        BAI_OUT="${BAM_OUT}.bai"
        if [ ! -f "${BAI_OUT}" ]; then
            echo "Indexing BAM for ${sample}..." >&2
            samtools index -@ ${THREADS} "${BAM_OUT}"
        fi
    fi
done

# --- 5. Variant Calling (lofreq) ---

for sample in "${SAMPLES[@]}"; do
    BAM_OUT="${RESULTS_DIR}/${sample}.bam"
    VCF_OUT_TEMP="${RESULTS_DIR}/${sample}.vcf"
    VCF_OUT_GZ="${RESULTS_DIR}/${sample}.vcf.gz"
    TBI_OUT="${RESULTS_DIR}/${sample}.vcf.gz.tbi"

    if [ ! -f "${BAM_OUT}" ] || [ ! -s "${BAM_OUT}" ]; then
        echo "Skipping variant calling for ${sample}: BAM file not found or empty." >&2
        continue
    fi

    # Check if VCF.gz and TBI exist (indicating successful prior run)
    if [ -f "${VCF_OUT_GZ}" ] && [ -f "${TBI_OUT}" ]; then
        echo "Skipping variant calling for ${sample}: VCF/TBI already exist." >&2
        continue
    fi

    echo "Calling variants for ${sample}..." >&2
    
    # 6. Variant calling with lofreq
    lofreq call-parallel --threads ${THREADS} \
        --reference "${REF_FA}" \
        --input "${BAM_OUT}" \
        --output "${VCF_OUT_TEMP}"

    if [ $? -ne 0 ]; then
        echo "Error during lofreq call-parallel for ${sample}." >&2
        exit 1
    fi

    # 7. VCF compression and indexing
    bgzip -c "${VCF_OUT_TEMP}" > "${VCF_OUT_GZ}"
    tabix -p vcf "${VCF_OUT_GZ}"
    
    # Cleanup intermediate file
    rm "${VCF_OUT_TEMP}"
done

# --- 8. Collapse Step ---

# Check if all VCFs exist before attempting collapse
ALL_VCFS_READY=true
for sample in "${SAMPLES[@]}"; do
    if [ ! -f "${RESULTS_DIR}/${sample}.vcf.gz" ]; then
        ALL_VCFS_READY=false
        break
    fi
done

if ! $ALL_VCFS_READY; then
    echo "Skipping collapse step: Not all per-sample VCFs are available." >&2
else
    # Check if the final TSV exists and is recent enough
    if [ -f "${COLLAPSED_TSV}" ] && [ $(find "${RESULTS_DIR}" -maxdepth 1 -type f -name "*.vcf.gz" -print0 | xargs -0 stat -c %Y | sort -n | tail -1) -ge $(stat -c %Y "${COLLAPSED_TSV}" 2>/dev/null || echo 0) ]; then
        echo "Skipping collapse step: ${COLLAPSED_TSV} is up to date." >&2
    else
        echo "Generating collapsed variant table..." >&2
        
        HEADER="sample\tchrom\tpos\tref\talt\taf"
        
        # Use printf to build the command string for bcftools query, ensuring the sample name is prepended
        # We use a temporary file to collect all sample outputs before piping to final concatenation
        TEMP_OUTPUT=$(mktemp)
        
        # Write header first
        echo -e "${HEADER}" > "${TEMP_OUTPUT}"
        
        # Append data for each sample
        for sample in "${SAMPLES[@]}"; do
            VCF_FILE="${RESULTS_DIR}/${sample}.vcf.gz"
            if [ -f "${VCF_FILE}" ]; then
                # bcftools query format: sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF
                bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "${VCF_FILE}" >> "${TEMP_OUTPUT}"
            fi
        done
        
        # Overwrite the final output
        mv "${TEMP_OUTPUT}" "${COLLAPSED_TSV}"
    fi
fi

exit 0