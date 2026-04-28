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
if [ ! -f "$REF_DIR/chrM.fa.fai" ]; then
    samtools faidx "$REF"
fi

if [ ! -f "$REF_DIR/chrM.fa.bwt" ]; then
    bwa index "$REF"
fi

# 3-7. Per-sample processing
for sample in "${SAMPLES[@]}"; do
    BAM="$RESULTS_DIR/${sample}.bam"
    BAI="$BAM.bai"
    VCF="$RESULTS_DIR/${sample}.vcf"
    VCF_GZ="${VCF}.gz"
    VCF_TBI="${VCF_GZ}.tbi"

    # Skip if all outputs exist and are newer than inputs
    if [ -f "$BAI" ] && [ -f "$VCF_TBI" ]; then
        # Check if inputs are newer
        need_work=0
        if [ "$RAW_DIR/${sample}_1.fq.gz" -nt "$BAI" ] || [ "$RAW_DIR/${sample}_2.fq.gz" -nt "$BAI" ] || [ "$REF" -nt "$BAI" ]; then
            need_work=1
        fi
        if [ "$need_work" -eq 0 ] && [ -f "$VCF_TBI" ]; then
            continue
        fi
    fi

    # 3. Alignment
    # Note: The prompt specifies exact -R format with literal \t.
    # We must ensure the shell does not interpret \t as a tab before passing to bwa.
    # Using single quotes for the RG string prevents shell interpretation.
    RG="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"

    # 4. SAM -> sorted BAM
    bwa mem -t "$THREADS" -R "$RG" "$REF" \
        "$RAW_DIR/${sample}_1.fq.gz" \
        "$RAW_DIR/${sample}_2.fq.gz" | \
        samtools sort -@ "$THREADS" -o "$BAM" -

    # 5. BAM indexing
    samtools index -@ "$THREADS" "$BAM"

    # 6. Variant calling with lofreq
    lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$VCF" "$BAM"

    # 7. VCF compression and indexing
    bgzip -f "$VCF"
    tabix -p vcf "$VCF_GZ"

    # Remove intermediate uncompressed VCF
    rm -f "$VCF"
done

# 8. Collapse step
COLLAPSED="$RESULTS_DIR/collapsed.tsv"
need_collapse=0

# Check if collapsed.tsv needs rebuilding
if [ -f "$COLLAPSED" ]; then
    for sample in "${SAMPLES[@]}"; do
        VCF_GZ="$RESULTS_DIR/${sample}.vcf.gz"
        if [ "$VCF_GZ" -nt "$COLLAPSED" ]; then
            need_collapse=1
            break
        fi
    done
else
    need_collapse=1
fi

if [ "$need_collapse" -eq 1 ]; then
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