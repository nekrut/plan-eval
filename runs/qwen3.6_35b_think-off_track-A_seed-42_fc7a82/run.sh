#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF="data/ref/chrM.fa"

mkdir -p results

# 2. Reference indexing (once)
if [ ! -f "data/ref/chrM.fa.fai" ] || [ ! -f "data/ref/chrM.fa.amb" ]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# 3-7. Per-sample processing
for sample in "${SAMPLES[@]}"; do
    BAM="results/${sample}.bam"
    BAI="results/${sample}.bam.bai"
    VCF_GZ="results/${sample}.vcf.gz"
    TBI="results/${sample}.vcf.gz.tbi"

    # Skip if BAM and BAI already exist
    if [ -f "$BAM" ] && [ -f "$BAI" ]; then
        :
    else
        # 3. Alignment
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            "$REF" \
            "data/raw/${sample}_1.fq.gz" \
            "data/raw/${sample}_2.fq.gz" | \
        # 4. Sort BAM
        samtools sort -@ "$THREADS" -o "$BAM" -

        # 5. Index BAM
        samtools index -@ "$THREADS" "$BAM"
    fi

    # Skip if VCF and TBI already exist
    if [ -f "$VCF_GZ" ] && [ -f "$TBI" ]; then
        :
    else
        # 6. Variant calling with lofreq
        lofreq call-parallel -f "$REF" -r "$BAM" -o "results/${sample}.vcf"

        # 7. Compress and index VCF
        bgzip -c "results/${sample}.vcf" > "$VCF_GZ"
        tabix -p vcf "$VCF_GZ"
        rm -f "results/${sample}.vcf"
    fi
done

# 8. Collapse step
COLLAPSED="results/collapsed.tsv"
REBUILD_COLLAPSED=0

# Check if collapsed.tsv needs rebuilding
if [ ! -f "$COLLAPSED" ]; then
    REBUILD_COLLAPSED=1
else
    for sample in "${SAMPLES[@]}"; do
        if [ "results/${sample}.vcf.gz" -nt "$COLLAPSED" ]; then
            REBUILD_COLLAPSED=1
            break
        fi
    done
fi

if [ "$REBUILD_COLLAPSED" -eq 1 ]; then
    # Write header
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED"
    
    # Append data for each sample
    for sample in "${SAMPLES[@]}"; do
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "results/${sample}.vcf.gz" >> "$COLLAPSED"
    done
fi