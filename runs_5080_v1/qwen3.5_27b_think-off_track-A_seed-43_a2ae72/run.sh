#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF="data/ref/chrM.fa"
REF_DIR="data/ref"
RAW_DIR="data/raw"
OUT_DIR="results"
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p "$OUT_DIR"

# 2. Reference indexing
if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.amb" ]]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# 3-7. Per-sample processing
for sample in "${SAMPLES[@]}"; do
    BAM="${OUT_DIR}/${sample}.bam"
    BAI="${OUT_DIR}/${sample}.bam.bai"
    VCF_GZ="${OUT_DIR}/${sample}.vcf.gz"
    TBI="${OUT_DIR}/${sample}.vcf.gz.tbi"
    R1="${RAW_DIR}/${sample}_1.fq.gz"
    R2="${RAW_DIR}/${sample}_2.fq.gz"

    # Idempotency check: if final index exists and is newer than BAM, skip
    if [[ -f "$TBI" ]] && [[ "$TBI" -nt "$BAM" ]] || [[ ! -f "$BAM" ]]; then
        # If BAM doesn't exist, we must run. If TBI exists but is older than BAM, we must run.
        # The logic below: if TBI exists AND is newer than BAM, skip. Otherwise run.
        if [[ -f "$TBI" ]] && [[ "$TBI" -nt "$BAM" ]]; then
            continue
        fi
    fi

    # 3 & 4. Alignment and sorting
    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" "$REF" "$R1" "$R2" | \
        samtools sort -@ "$THREADS" -o "$BAM" -

    # 5. BAM indexing
    samtools index -@ "$THREADS" "$BAM"

    # 6. Variant calling
    VCF_TMP="${OUT_DIR}/${sample}.vcf"
    lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$VCF_TMP" "$BAM"

    # 7. VCF compression and indexing
    bgzip -c "$VCF_TMP" > "$VCF_GZ"
    rm -f "$VCF_TMP"
    tabix -p vcf "$VCF_GZ"
done

# 8. Collapse step
COLLAPSED="${OUT_DIR}/collapsed.tsv"
# Check if any VCF is newer than the collapsed file
REBUILD=false
if [[ ! -f "$COLLAPSED" ]]; then
    REBUILD=true
else
    for sample in "${SAMPLES[@]}"; do
        VCF_GZ="${OUT_DIR}/${sample}.vcf.gz"
        if [[ "$VCF_GZ" -nt "$COLLAPSED" ]]; then
            REBUILD=true
            break
        fi
    done
fi

if [[ "$REBUILD" == true ]]; then
    {
        printf "sample\tchrom\tpos\tref\talt\taf\n"
        for sample in "${SAMPLES[@]}"; do
            VCF_GZ="${OUT_DIR}/${sample}.vcf.gz"
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ"
        done
    } > "$COLLAPSED"
fi