#!/usr/bin/env bash
set -euo pipefail

# Ensure reference is indexed
if [ ! -f data/ref/chrM.fa.fai ]; then
    samtools faidx data/ref/chrM.fa
fi

samples=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

for sample in "${samples[@]}"; do
    bam="results/${sample}.bam"
    bai="results/${sample}.bam.bai"
    vcf="results/${sample}.vcf.gz"
    tbi="results/${sample}.vcf.gz.tbi"
    filtered_vcf="results/${sample}_filtered.vcf"
    annotated_vcf="results/${sample}_annotated.vcf"
    temp_file="results/${sample}_temp.tsv"

    # Skip if BAM and BAI exist
    if [ -f "$bam" ] && [ -f "$bai" ]; then
        continue
    fi

    # Align with bwa
    if [ ! -f "$bam" ]; then
        bwa mem -t 4 data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
        samtools sort -o "$bam" -@ 4
        samtools index "$bam" "$bai"
    fi

    # Generate VCF
    if [ ! -f "$vcf" ]; then
        bcftools mpileup -f data/ref/chrM.fa "$bam" | \
        bcftools call -mv -Ov -o "$filtered_vcf" -@ 4
        bcftools filter -i 'TYPE != snp' "$filtered_vcf" -o "$annotated_vcf" -@ 4
        bcftools view -Oz -o "$vcf" "$annotated_vcf"
        tabix -p vcf "$vcf"
    fi
done

# Process VCFs to generate collapsed.tsv
collapsed_tsv="results/collapsed.tsv"
temp_dir="results/temp"
mkdir -p "$temp_dir"

for sample in "${samples[@]}"; do
    vcf="results/${sample}.vcf.gz"
    temp_file="$temp_dir/${sample}_temp.tsv"

    if [ -f "$temp_file" ]; then
        continue
    fi

    # Extract sample name
    sample_name="$sample"

    # Add AF with SnpSift and extract fields
    snpSift expr -e 'AF = (FORMAT/GT == "0/1" || FORMAT/GT == "1/0" ? 0.5 : (FORMAT/GT == "1/1" ? 1 : 0))' "$vcf" | \
    awk -v sample="$sample_name" 'BEGIN { OFS="\t"; print "sample", "chrom", "pos", "ref", "alt", "af"} { if ($0 ~ /^#/) next; chrom = $1; pos = $2; ref = $3; alt = $4; af = $8; print sample, chrom, pos, ref, alt, af }' > "$temp_file"
done

# Combine all temp files into collapsed.tsv
cat "$temp_dir"/*.tsv | sort -k2,2n -k3,3n > "$collapsed_tsv"