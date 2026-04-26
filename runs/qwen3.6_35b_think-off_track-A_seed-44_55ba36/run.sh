#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF="data/ref/chrM.fa"
REF_DIR="data/ref"
RAW_DIR="data/raw"
RES_DIR="results"

mkdir -p "$RES_DIR"

# 2. Reference indexing (once)
if [ ! -f "$REF_DIR/chrM.fa.fai" ]; then
    samtools faidx "$REF"
fi
if [ ! -f "$REF_DIR/chrM.fa.bwt" ]; then
    bwa index "$REF"
fi

# 3-7. Per-sample processing
for sample in "${SAMPLES[@]}"; do
    BAM="$RES_DIR/${sample}.bam"
    BAI="$BAM.bai"
    VCF="$RES_DIR/${sample}.vcf"
    VCF_GZ="$VCF.gz"
    TBI="$VCF_GZ.tbi"

    # 3-4. Alignment and sorting
    if [ ! -f "$BAM" ] || [ "$RAW_DIR/${sample}_1.fq.gz" -nt "$BAM" ] || [ "$RAW_DIR/${sample}_2.fq.gz" -nt "$BAM" ]; then
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            "$REF" \
            "$RAW_DIR/${sample}_1.fq.gz" \
            "$RAW_DIR/${sample}_2.fq.gz" \
        | samtools sort -@ "$THREADS" -o "$BAM" -
    fi

    # 5. BAM indexing
    if [ ! -f "$BAI" ] || [ "$BAM" -nt "$BAI" ]; then
        samtools index -@ "$THREADS" "$BAM"
    fi

    # 6. Variant calling with lofreq
    if [ ! -f "$VCF" ] || [ "$BAM" -nt "$VCF" ]; then
        lofreq call-parallel -f "$REF" -d -o "$VCF" -r 1-16569 "$BAM"
    fi

    # 7. VCF compression and indexing
    if [ ! -f "$VCF_GZ" ] || [ "$VCF" -nt "$VCF_GZ" ]; then
        bgzip -c "$VCF" > "$VCF_GZ"
        tabix -p vcf "$VCF_GZ"
        rm -f "$VCF"
    fi
done

# 8. Collapse step
COLLAPSED="$RES_DIR/collapsed.tsv"
HEADER="sample	chrom	pos	ref	alt	af"

# Check if we need to regenerate
NEED_REBUILD=false
if [ ! -f "$COLLAPSED" ]; then
    NEED_REBUILD=true
else
    for sample in "${SAMPLES[@]}"; do
        if [ "$RES_DIR/${sample}.vcf.gz" -nt "$COLLAPSED" ]; then
            NEED_REBUILD=true
            break
        fi
    done
fi

if [ "$NEED_REBUILD" = true ]; then
    {
        echo -e "$HEADER"
        for sample in "${SAMPLES[@]}"; do
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$RES_DIR/${sample}.vcf.gz"
        done
    } > "$COLLAPSED"
fi