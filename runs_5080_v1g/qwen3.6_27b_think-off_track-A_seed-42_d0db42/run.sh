#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF="data/ref/chrM.fa"
RAW_DIR="data/raw"
OUT_DIR="results"

mkdir -p "${OUT_DIR}"

# Step 2: Reference indexing (idempotent)
if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.bwt" ]]; then
    samtools faidx "${REF}"
    bwa index "${REF}"
fi

# Step 3-7: Per-sample alignment, sorting, indexing, calling, compression
for sample in "${SAMPLES[@]}"; do
    bam_file="${OUT_DIR}/${sample}.bam"
    bai_file="${OUT_DIR}/${sample}.bam.bai"
    vcf_gz_file="${OUT_DIR}/${sample}.vcf.gz"
    tbi_file="${OUT_DIR}/${sample}.vcf.gz.tbi"

    # Check if final VCF index exists; if so, skip all steps for this sample
    if [[ -f "${tbi_file}" ]]; then
        continue
    fi

    # Step 3 & 4: Alignment and sorting
    # bwa mem with read group, piped to samtools sort
    # Read group string uses literal \t as required by bwa mem
    rg_str="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
    
    bwa mem -t "${THREADS}" -R "${rg_str}" "${REF}" \
        "${RAW_DIR}/${sample}_1.fq.gz" \
        "${RAW_DIR}/${sample}_2.fq.gz" | \
        samtools sort -@ "${THREADS}" -o "${bam_file}" -

    # Step 5: BAM indexing
    samtools index -@ "${THREADS}" "${bam_file}"

    # Step 6: Variant calling with lofreq call-parallel
    # Note: lofreq call-parallel outputs to stdout by default if --out is not specified, 
    # but the plan specifies --out results/{sample}.vcf. We follow the plan.
    lofreq call-parallel --pp-threads "${THREADS}" \
        --ref "${REF}" \
        --out "${OUT_DIR}/${sample}.vcf" \
        "${bam_file}"

    # Step 7: VCF compression and indexing
    # Compress with bgzip
    bgzip -c "${OUT_DIR}/${sample}.vcf" > "${vcf_gz_file}"
    
    # Index with tabix
    tabix -p vcf "${vcf_gz_file}"
    
    # Remove intermediate uncompressed VCF
    rm -f "${OUT_DIR}/${sample}.vcf"
done

# Step 8: Collapse step -> results/collapsed.tsv
# Rebuild only if any input VCF is newer than the TSV or if TSV doesn't exist
tsv_file="${OUT_DIR}/collapsed.tsv"
needs_rebuild=false

if [[ ! -f "${tsv_file}" ]]; then
    needs_rebuild=true
else
    for sample in "${SAMPLES[@]}"; do
        vcf_gz="${OUT_DIR}/${sample}.vcf.gz"
        if [[ "${vcf_gz}" -nt "${tsv_file}" ]]; then
            needs_rebuild=true
            break
        fi
    done
fi

if [[ "${needs_rebuild}" == true ]]; then
    # Create header
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "${tsv_file}"
    
    # Append data from each sample
    for sample in "${SAMPLES[@]}"; do
        vcf_gz="${OUT_DIR}/${sample}.vcf.gz"
        # bcftools query format: sample name (literal), chrom, pos, ref, alt, AF from INFO
        # Note: %INFO/AF might be empty if AF is not present, but lofreq usually provides it.
        # If AF is missing, bcftools might output empty string.
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "${vcf_gz}" >> "${tsv_file}"
    done
fi