#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF="data/ref/chrM.fa"
REF_DIR="data/ref"
RAW_DIR="data/raw"
RESULTS_DIR="results"

mkdir -p "$RESULTS_DIR"

# 2. Reference indexing (once)
if [ ! -f "$REF_DIR/chrM.fa.fai" ] || [ ! -f "$REF_DIR/chrM.fa.bwt" ]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# 3-7. Per-sample processing
for sample in "${SAMPLES[@]}"; do
    BAM="$RESULTS_DIR/${sample}.bam"
    BAI="$BAM.bai"
    VCF="$RESULTS_DIR/${sample}.vcf"
    VCF_GZ="$VCF.gz"
    VCF_TBI="$VCF_GZ.tbi"

    # Skip if all outputs exist
    if [ -f "$BAI" ] && [ -f "$VCF_TBI" ]; then
        continue
    fi

    # 3. Alignment
    # Note: \t is passed as literal backslash-t for bwa to parse
    bwa mem -t "$THREADS" \
        -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
        "$REF" \
        "${RAW_DIR}/${sample}_1.fq.gz" \
        "${RAW_DIR}/${sample}_2.fq.gz" \
        | samtools sort -@ "$THREADS" -o "$BAM" -

    # 5. BAM indexing
    samtools index -@ "$THREADS" "$BAM"

    # 6. Variant calling with lofreq
    lofreq call-parallel \
        -f "$REF" \
        -r "$BAM" \
        -o "$VCF" \
        --pp-threads "$THREADS" \
        --no-indels

    # 7. VCF compression and indexing
    bgzip -f "$VCF"
    tabix -p vcf "$VCF_GZ"
done

# 8. Collapse step
COLLAPSED="$RESULTS_DIR/collapsed.tsv"
REBUILD_COLLAPSED=0

# Check if any VCF_GZ is newer than COLLAPSED or if COLLAPSED doesn't exist
for sample in "${SAMPLES[@]}"; do
    VCF_GZ="$RESULTS_DIR/${sample}.vcf.gz"
    if [ ! -f "$COLLAPSED" ]; then
        REBUILD_COLLAPSED=1
        break
    fi
    if [ "$VCF_GZ" -nt "$COLLAPSED" ]; then
        REBUILD_COLLAPSED=1
        break
    fi
done

if [ "$REBUILD_COLLAPSED" -eq 1 ]; then
    # Write header
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED"
    
    # Append data for each sample
    for sample in "${SAMPLES[@]}"; do
        VCF_GZ="$RESULTS_DIR/${sample}.vcf.gz"
        if [ -f "$VCF_GZ" ]; then
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ" >> "$COLLAPSED"
        fi
    done
fi