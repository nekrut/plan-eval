#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF_DIR="data/ref"
RAW_DIR="data/raw"
RESULTS_DIR="results"
REF_FASTA="$REF_DIR/chrM.fa"
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p "$RESULTS_DIR"

# Index reference genome
if [ ! -f "$REF_FASTA.bwt" ]; then
    bwa index "$REF_FASTA"
fi

if [ ! -f "$REF_FASTA.fai" ]; then
    samtools faidx "$REF_FASTA"
fi

# Process each sample
for sample in "${SAMPLES[@]}"; do
    R1="$RAW_DIR/${sample}_1.fq.gz"
    R2="$RAW_DIR/${sample}_2.fq.gz"
    BAM="$RESULTS_DIR/${sample}.bam"
    BAI="$RESULTS_DIR/${sample}.bam.bai"
    VCF_GZ="$RESULTS_DIR/${sample}.vcf.gz"
    VCF_TBI="$RESULTS_DIR/${sample}.vcf.gz.tbi"
    
    # Skip if all outputs exist
    if [ -f "$BAM" ] && [ -f "$BAI" ] && [ -f "$VCF_GZ" ] && [ -f "$VCF_TBI" ]; then
        continue
    fi
    
    # Map paired-end reads
    SAM="$RESULTS_DIR/${sample}.sam"
    if [ ! -f "$BAM" ]; then
        bwa mem -t "$THREADS" "$REF_FASTA" "$R1" "$R2" > "$SAM"
        samtools view -b -S "$SAM" | samtools sort -@ "$THREADS" -o "$BAM"
        rm "$SAM"
    fi
    
    # Index BAM
    if [ ! -f "$BAI" ]; then
        samtools index "$BAM"
    fi
    
    # Call variants
    VCF="$RESULTS_DIR/${sample}.vcf"
    if [ ! -f "$VCF_GZ" ]; then
        lofreq call -f "$REF_FASTA" -o "$VCF" "$BAM"
        bgzip -f "$VCF"
    fi
    
    # Index VCF
    if [ ! -f "$VCF_TBI" ]; then
        tabix -p vcf "$VCF_GZ"
    fi
done

# Create collapsed TSV
COLLAPSED="$RESULTS_DIR/collapsed.tsv"
if [ ! -f "$COLLAPSED" ]; then
    {
        echo -e "sample\tchrom\tpos\tref\talt\taf"
        for sample in "${SAMPLES[@]}"; do
            VCF_GZ="$RESULTS_DIR/${sample}.vcf.gz"
            if [ -f "$VCF_GZ" ]; then
                zcat "$VCF_GZ" | grep -v '^#' | while IFS=$'\t' read -r chrom pos id ref alt qual filter info rest; do
                    af="0"
                    if [[ "$info" =~ AF=([0-9.]+) ]]; then
                        af="${BASH_REMATCH[1]}"
                    fi
                    echo -e "$sample\t$chrom\t$pos\t$ref\t$alt\t$af"
                done
            fi
        done
    } > "$COLLAPSED"
fi