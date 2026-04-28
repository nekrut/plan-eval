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

# Check if all required outputs exist to determine if work is needed
ALL_OUTPUTS_EXIST=true
for sample in "${SAMPLES[@]}"; do
    if [[ ! -f "${RESULTS_DIR}/${sample}.vcf.gz.tbi" ]]; then
        ALL_OUTPUTS_EXIST=false
        break
    fi
done

if $ALL_OUTPUTS_EXIST && [[ -f "${COLLAPSED_FILE}" ]]; then
    # Check if the collapsed file is up-to-date relative to the VCFs
    latest_vcf_time=0
    for sample in "${SAMPLES[@]}"; do
        vcf_file="${RESULTS_DIR}/${sample}.vcf.gz"
        if [[ -f "$vcf_file" ]]; then
            vcf_time=$(stat -c %Y "$vcf_file")
            if (( vcf_time > latest_vcf_time )); then
                latest_vcf_time=$vcf_time
            fi
        fi
    done

    if [[ -f "${COLLAPSED_FILE}" ]]; then
        tsv_time=$(stat -c %Y "${COLLAPSED_FILE}")
        if (( tsv_time >= latest_vcf_time )); then
            : # Outputs are up to date
        else
            echo "Warning: Rebuilding ${COLLAPSED_FILE} because it is older than the input VCFs." >&2
            # Proceed to step 8 if necessary
        fi
    else
        echo "Warning: Rebuilding ${COLLAPSED_FILE} because it does not exist." >&2
        # Proceed to step 8
    fi
fi

# --- 2. Reference Indexing ---
echo "Indexing reference genome..." >&2
if [[ ! -f "${REF_FA}.fai" ]]; then
    samtools faidx "${REF_FA}"
fi

if [[ ! -f "${REF_FA}.bwt" ]]; then
    bwa index "${REF_FA}"
fi

# --- 3. Per-sample alignment with bwa mem ---
for sample in "${SAMPLES[@]}"; do
    echo "Processing sample: ${sample}" >&2
    R1="${RAW_DIR}/${sample}_1.fq.gz"
    R2="${RAW_DIR}/${sample}_2.fq.gz"
    OUTPUT_BAM="${RESULTS_DIR}/${sample}.bam"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "Skipping ${sample}: Input files not found." >&2
        continue
    fi

    # Check if BAM exists and is recent enough to skip alignment
    if [[ -f "$OUTPUT_BAM" ]] && [[ -f "${OUTPUT_BAM}.bai" ]]; then
        # Simple check: if BAM exists, assume it's done unless we explicitly need to re-run
        # For strict idempotency, we rely on the subsequent steps failing if the BAM is missing.
        :
    fi

    # Read Group construction: -R "@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA"
    # Note: The literal backslash-t must be passed carefully.
    RG_LINE="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"

    # Run bwa mem and pipe to samtools sort
    if ! bwa mem -t ${THREADS} -R "${RG_LINE}" "${REF_FA}" "${R1}" "${R2}" 2> /dev/null | \
        samtools sort -@ ${THREADS} -o "${OUTPUT_BAM}" -T "${RESULTS_DIR}/${sample}.tmp.sorted.bam"; then
        echo "Error during bwa mem/samtools sort for ${sample}. Skipping." >&2
        rm -f "${RESULTS_DIR}/${sample}.tmp.sorted.bam"
        continue
    fi
    mv "${RESULTS_DIR}/${sample}.tmp.sorted.bam" "${OUTPUT_BAM}"
done

# --- 4. BAM Indexing ---
for sample in "${SAMPLES[@]}"; do
    BAM_FILE="${RESULTS_DIR}/${sample}.bam"
    BAI_FILE="${RESULTS_DIR}/${sample}.bam.bai"
    if [[ -f "$BAM_FILE" ]] && [[ ! -f "$BAI_FILE" ]]; then
        echo "Indexing BAM for ${sample}..." >&2
        samtools index -@ ${THREADS} "${BAM_FILE}"
    fi
done

# --- 5. Variant calling with lofreq call-parallel ---
for sample in "${SAMPLES[@]}"; do
    BAM_FILE="${RESULTS_DIR}/${sample}.bam"
    VCF_GZ="${RESULTS_DIR}/${sample}.vcf.gz"
    VCF_TBI="${RESULTS_DIR}/${sample}.vcf.gz.tbi"

    if [[ ! -f "$BAM_FILE" ]]; then
        echo "Skipping variant calling for ${sample}: BAM file missing." >&2
        continue
    fi

    # Check if VCF.tbi exists, implying VCF.gz exists
    if [[ -f "$VCF_TBI" ]]; then
        echo "Skipping variant calling for ${sample}: VCF index exists." >&2
        continue
    fi

    echo "Calling variants for ${sample}..." >&2
    # lofreq call-parallel --pp-threads 4 --verbose --ref data/ref/chrM.fa --out results/{sample}.vcf --sig --bonf results/{sample}.bam
    if ! lofreq call-parallel --pp-threads ${THREADS} \
        --ref "${REF_FA}" \
        --out "${RESULTS_DIR}/${sample}.vcf" \
        --sig \
        --bonf \
        "${BAM_FILE}"; then
        echo "Error during lofreq call-parallel for ${sample}. Skipping." >&2
        continue
    fi

    # --- 6. VCF compression and indexing ---
    VCF_TEMP="${RESULTS_DIR}/${sample}.vcf"
    VCF_GZ="${RESULTS_DIR}/${sample}.vcf.gz"

    if [[ ! -f "$VCF_TEMP" ]]; then
        echo "Error: Intermediate VCF file not found for ${sample}." >&2
        continue
    fi

    if [[ ! -f "$VCF_GZ" ]]; then
        echo "Compressing VCF for ${sample}..." >&2
        bgzip -c "${VCF_TEMP}" > "${VCF_GZ}"
    fi

    if [[ ! -f "$VCF_TBI" ]]; then
        echo "Indexing VCF for ${sample}..." >&2
        tabix -p vcf "${VCF_GZ}"
    fi

    # Cleanup intermediate file
    rm -f "${VCF_TEMP}"
done

# --- 7. Collapse step → results/collapsed.tsv ---
echo "Collapsing results into ${COLLAPSED_FILE}..." >&2

# Check if any VCF exists to proceed
if ! ls "${RESULTS_DIR}"/*.vcf.gz 1> /dev/null 2>&1; then
    echo "No VCF files found to collapse." >&2
    exit 0
fi

# Determine if rebuilding is necessary based on timestamp check at the start
if [[ ! -f "$COLLAPSED_FILE" ]] || [[ $(stat -c %Y "${COLLAPSED_FILE}") -lt $(stat -c %Y "${RESULTS_DIR}/M117-bl.vcf.gz") ]]; then
    
    HEADER="sample\tchrom\tpos\tref\talt\taf"
    
    # Use printf to handle the header and then loop for data
    {
        echo -e "${HEADER}"
        for sample in "${SAMPLES[@]}"; do
            VCF_FILE="${RESULTS_DIR}/${sample}.vcf.gz"
            if [[ -f "$VCF_FILE" ]]; then
                # bcftools query format: {sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n
                bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "${VCF_FILE}"
            fi
        done
    } | sort -k1,1 -k2,2 -k3,3 > "${COLLAPSED_FILE}"
fi

exit 0