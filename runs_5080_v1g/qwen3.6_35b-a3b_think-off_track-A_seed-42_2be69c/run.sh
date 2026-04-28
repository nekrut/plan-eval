#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF="data/ref/chrM.fa"
REF_IDX="data/ref/chrM.fa.fai"
REF_BWA="data/ref/chrM.fa"

mkdir -p results

# Step 2: Reference indexing (idempotent)
if [ ! -f "$REF_IDX" ]; then
    samtools faidx "$REF"
fi
if [ ! -f "${REF_BWA}.amb" ]; then
    bwa index "$REF"
fi

# Step 3-7: Per-sample processing
for sample in "${SAMPLES[@]}"; do
    BAM="results/${sample}.bam"
    BAI="results/${sample}.bam.bai"
    VCF="results/${sample}.vcf"
    VCF_GZ="results/${sample}.vcf.gz"
    VCF_TBI="results/${sample}.vcf.gz.tbi"

    # Step 3-4: Alignment and sorting
    if [ ! -f "$BAM" ] || [ "$REF" -nt "$BAM" ] || [ "data/raw/${sample}_1.fq.gz" -nt "$BAM" ] || [ "data/raw/${sample}_2.fq.gz" -nt "$BAM" ]; then
        bwa mem -t $THREADS \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            "$REF" \
            "data/raw/${sample}_1.fq.gz" \
            "data/raw/${sample}_2.fq.gz" \
        | samtools sort -@ $THREADS -o "$BAM" -
    fi

    # Step 5: BAM indexing
    if [ ! -f "$BAI" ]; then
        samtools index -@ $THREADS "$BAM"
    fi

    # Step 6: Variant calling
    if [ ! -f "$VCF" ] || [ "$BAM" -nt "$VCF" ]; then
        lofreq call-parallel --pp-threads $THREADS --verbose \
            --ref "$REF" \
            --out "$VCF" \
            --sig \
            --bonf \
            "$BAM"
    fi

    # Step 7: VCF compression and indexing
    if [ ! -f "$VCF_GZ" ] || [ "$VCF" -nt "$VCF_GZ" ]; then
        bgzip -f "$VCF"
        tabix -p vcf "$VCF_GZ"
    fi
done

# Step 8: Collapse step
COLLAPSED="results/collapsed.tsv"
HEADER="sample	chrom	pos	ref	alt	af"

# Check if collapsed.tsv needs rebuilding
NEED_REBUILD=false
if [ ! -f "$COLLAPSED" ]; then
    NEED_REBUILD=true
else
    for sample in "${SAMPLES[@]}"; do
        VCF_GZ="results/${sample}.vcf.gz"
        if [ "$VCF_GZ" -nt "$COLLAPSED" ]; then
            NEED_REBUILD=true
            break
        fi
    done
fi

if [ "$NEED_REBUILD" = true ]; then
    {
        echo -e "$HEADER"
        for sample in "${SAMPLES[@]}"; do
            VCF_GZ="results/${sample}.vcf.gz"
            if [ -f "$VCF_GZ" ]; then
                bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ"
            fi
        done
    } > "$COLLAPSED"
fi