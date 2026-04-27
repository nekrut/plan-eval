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

REF_INDEXES=( "$REF_DIR/chrM.fa.amb" "$REF_DIR/chrM.fa.ann" "$REF_DIR/chrM.fa.bwt" "$REF_DIR/chrM.fa.pac" "$REF_DIR/chrM.fa.sa" )
need_bwa=0
for idx in "${REF_INDEXES[@]}"; do
    if [ ! -f "$idx" ]; then
        need_bwa=1
        break
    fi
done

if [ "$need_bwa" -eq 1 ]; then
    bwa index "$REF"
fi

# 3-7. Per-sample processing
for sample in "${SAMPLES[@]}"; do
    bam="$RESULTS_DIR/${sample}.bam"
    bai="$bam.bai"
    vcf_gz="$RESULTS_DIR/${sample}.vcf.gz"
    vcf_tbi="$vcf_gz.tbi"

    # Skip if all outputs exist and are newer than inputs
    if [ -f "$vcf_tbi" ]; then
        # Check if any input is newer
        inputs_updated=0
        if [ "$RAW_DIR/${sample}_1.fq.gz" -nt "$vcf_tbi" ] || [ "$RAW_DIR/${sample}_2.fq.gz" -nt "$vcf_tbi" ]; then
            inputs_updated=1
        fi
        if [ "$REF" -nt "$vcf_tbi" ]; then
            inputs_updated=1
        fi
        if [ "$REF_DIR/chrM.fa.fai" -nt "$vcf_tbi" ]; then
            inputs_updated=1
        fi
        for idx in "${REF_INDEXES[@]}"; do
            if [ "$idx" -nt "$vcf_tbi" ]; then
                inputs_updated=1
                break
            fi
        done

        if [ "$inputs_updated" -eq 0 ]; then
            continue
        fi
    fi

    # 3. Alignment
    samtools fastq -1 /dev/null -2 /dev/null -s /dev/null -n "$bam" 2>/dev/null || true # dummy check if bam exists, but we need to generate it
    
    # Run bwa mem and sort
    bwa mem -t "$THREADS" \
        -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
        "$REF" \
        "$RAW_DIR/${sample}_1.fq.gz" \
        "$RAW_DIR/${sample}_2.fq.gz" | \
        samtools sort -@ "$THREADS" -o "$bam" -

    # 5. BAM indexing
    samtools index -@ "$THREADS" "$bam"

    # 6. Variant calling with lofreq
    vcf_plain="$RESULTS_DIR/${sample}.vcf"
    lofreq call-parallel -f "$REF" -i "$bam" -o "$vcf_plain" --pp-threads "$THREADS" --no-indels

    # 7. VCF compression and indexing
    bgzip -f "$vcf_plain"
    tabix -p vcf "$vcf_gz"

    # Cleanup intermediate
    rm -f "$vcf_plain"
done

# 8. Collapse step
collapsed="$RESULTS_DIR/collapsed.tsv"
need_collapse=0

# Check if collapsed.tsv needs rebuilding
if [ -f "$collapsed" ]; then
    for sample in "${SAMPLES[@]}"; do
        vcf_gz="$RESULTS_DIR/${sample}.vcf.gz"
        if [ "$vcf_gz" -nt "$collapsed" ]; then
            need_collapse=1
            break
        fi
    done
else
    need_collapse=1
fi

if [ "$need_collapse" -eq 1 ]; then
    # Write header
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$collapsed"
    
    # Append data
    for sample in "${SAMPLES[@]}"; do
        vcf_gz="$RESULTS_DIR/${sample}.vcf.gz"
        if [ -f "$vcf_gz" ]; then
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$vcf_gz" >> "$collapsed"
        fi
    done
fi