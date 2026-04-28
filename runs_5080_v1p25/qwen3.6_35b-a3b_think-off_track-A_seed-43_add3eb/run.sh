#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF="data/ref/chrM.fa"
REF_IDX="data/ref/chrM.fa.fai"
REF_BWA="data/ref/chrM.fa.amb"

mkdir -p results

# 2. Reference indexing (once)
if [ ! -f "$REF_IDX" ] || [ ! -f "$REF_BWA" ]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# 3-7. Per-sample processing
for sample in "${SAMPLES[@]}"; do
    BAM="results/${sample}.bam"
    BAI="results/${sample}.bam.bai"
    VCF="results/${sample}.vcf"
    VCF_GZ="results/${sample}.vcf.gz"
    VCF_TBI="results/${sample}.vcf.gz.tbi"

    # Skip if BAM and BAI already exist (alignment done)
    if [ ! -f "$BAI" ]; then
        # 3. Alignment
        bwa mem -t $THREADS \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            "$REF" \
            "data/raw/${sample}_1.fq.gz" \
            "data/raw/${sample}_2.fq.gz" | \
        # 4. Sort BAM
        samtools sort -@ $THREADS -o "$BAM" -

        # 5. Index BAM
        samtools index -@ $THREADS "$BAM"
    fi

    # Skip if VCF.GZ and TBI already exist (variant calling done)
    if [ ! -f "$VCF_TBI" ]; then
        # 6. Variant calling
        lofreq call-parallel --pp-threads $THREADS -f "$REF" -o "$VCF" "$BAM"

        # 7. Compress and index VCF
        bgzip -f "$VCF"
        tabix -p vcf "$VCF_GZ"
    fi
done

# 8. Collapse step
COLLAPSED="results/collapsed.tsv"
HEADER="sample	chrom	pos	ref	alt	af"

# Check if collapsed.tsv exists and is newer than all input VCF.GZ files
needs_collapse=0
if [ ! -f "$COLLAPSED" ]; then
    needs_collapse=1
else
    for sample in "${SAMPLES[@]}"; do
        if [ "results/${sample}.vcf.gz" -nt "$COLLAPSED" ]; then
            needs_collapse=1
            break
        fi
    done
fi

if [ "$needs_collapse" -eq 1 ]; then
    # Write header
    printf '%s\n' "$HEADER" > "$COLLAPSED"
    
    # Append data from each sample
    for sample in "${SAMPLES[@]}"; do
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "results/${sample}.vcf.gz" >> "$COLLAPSED"
    done
fi