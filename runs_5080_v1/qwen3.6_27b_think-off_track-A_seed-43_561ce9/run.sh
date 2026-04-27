#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF="data/ref/chrM.fa"
RAW_DIR="data/raw"
RES_DIR="results"

mkdir -p "$RES_DIR"

# 2. Reference indexing (idempotent)
if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.bwt" ]]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# 3-7. Per-sample alignment, sorting, indexing, calling, compression
for sample in "${SAMPLES[@]}"; do
    bam_file="${RES_DIR}/${sample}.bam"
    bai_file="${RES_DIR}/${sample}.bam.bai"
    vcf_gz_file="${RES_DIR}/${sample}.vcf.gz"
    tbi_file="${RES_DIR}/${sample}.vcf.gz.tbi"

    # Check if all downstream artifacts exist and are up-to-date
    # We check the final TBI file as the sentinel for the whole pipeline for this sample
    if [[ -f "$tbi_file" ]]; then
        # Verify dependencies haven't changed (simple timestamp check)
        # If BAM is newer than TBI, or VCF is newer than TBI, we might need to rebuild.
        # However, for strict idempotency without redoing work, we assume if TBI exists, we are good.
        # To be safer, we check if the BAM exists. If BAM is missing but TBI is, something is wrong, but we skip.
        # If BAM is newer than BAI, rebuild BAI.
        if [[ "$bam_file" -nt "$bai_file" ]]; then
            samtools index -@ "$THREADS" "$bam_file"
        fi
        # If VCF.gz is newer than TBI, rebuild TBI
        if [[ "$vcf_gz_file" -nt "$tbi_file" ]]; then
            tabix -p vcf "$vcf_gz_file"
        fi
        continue
    fi

    # 3. Alignment
    # 4. Sorting
    # Note: bwa mem stdout -> samtools sort
    # Read group string: literal \t
    rg="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
    
    bwa mem -t "$THREADS" -R "$rg" "$REF" \
        "${RAW_DIR}/${sample}_1.fq.gz" \
        "${RAW_DIR}/${sample}_2.fq.gz" | \
        samtools sort -@ "$THREADS" -o "$bam_file"

    # 5. BAM indexing
    samtools index -@ "$THREADS" "$bam_file"

    # 6. Variant calling with lofreq call-parallel
    # Output to uncompressed VCF first
    vcf_tmp="${RES_DIR}/${sample}.vcf"
    lofreq call-parallel \
        --ref "$REF" \
        --pp-threads "$THREADS" \
        -o "$vcf_tmp" \
        "$bam_file"

    # 7. VCF compression and indexing
    bgzip -c "$vcf_tmp" > "$vcf_gz_file"
    tabix -p vcf "$vcf_gz_file"
    
    # Remove intermediate uncompressed VCF
    rm -f "$vcf_tmp"
done

# 8. Collapse step -> results/collapsed.tsv
# Rebuild only if any input VCF is newer than the TSV
NEED_COLLAPSE=false
for sample in "${SAMPLES[@]}"; do
    vcf_gz_file="${RES_DIR}/${sample}.vcf.gz"
    if [[ ! -f "$RES_DIR/collapsed.tsv" ]] || [[ "$vcf_gz_file" -nt "$RES_DIR/collapsed.tsv" ]]; then
        NEED_COLLAPSE=true
        break
    fi
done

if [[ "$NEED_COLLAPSE" == true ]]; then
    # Create header
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$RES_DIR/collapsed.tsv"
    
    # Append data from each sample
    for sample in "${SAMPLES[@]}"; do
        vcf_gz_file="${RES_DIR}/${sample}.vcf.gz"
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$vcf_gz_file" >> "$RES_DIR/collapsed.tsv"
    done
fi