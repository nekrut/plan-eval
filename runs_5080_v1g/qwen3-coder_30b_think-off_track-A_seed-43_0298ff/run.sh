#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

mkdir -p results
cd "$(dirname "$0")"

REF=data/ref/chrM.fa
REF_FAIDX=${REF}.fai
REF_BWA_INDEXES=(${REF}.amb ${REF}.ann ${REF}.bwt ${REF}.pac ${REF}.sa)

# Index reference if needed
if [[ ! -f "$REF_FAIDX" ]] || [[ ! -all ${REF_BWA_INDEXES[@]} ]]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# Process each sample
for sample in "${SAMPLES[@]}"; do
    BAM=results/${sample}.bam
    BAM_BAI=results/${sample}.bam.bai
    VCF=results/${sample}.vcf
    VCF_GZ=results/${sample}.vcf.gz
    VCF_GZ_TBI=results/${sample}.vcf.gz.tbi

    # Skip if all outputs exist and are newer than inputs
    if [[ -f "$VCF_GZ_TBI" ]] && [[ "$VCF_GZ_TBI" -nt "$BAM" ]] && [[ "$VCF_GZ_TBI" -nt "${REF}" ]]; then
        continue
    fi

    # Align with bwa mem
    bwa mem -t "$THREADS" \
        -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
        "$REF" \
        data/raw/${sample}_1.fq.gz \
        data/raw/${sample}_2.fq.gz \
        | samtools sort -@ "$THREADS" -o "$BAM"

    # Index BAM
    samtools index -@ "$THREADS" "$BAM"

    # Call variants with lofreq
    lofreq call-parallel --pp-threads "$THREADS" --verbose \
        --ref "$REF" --out "$VCF" \
        --sig \
        --bonf \
        "$BAM"

    # Compress and index VCF
    bgzip "$VCF"
    tabix -p vcf "$VCF_GZ"

    # Remove uncompressed VCF
    rm -f "$VCF"
done

# Collapse all VCFs into single TSV
COLLAPSED=results/collapsed.tsv
if [[ ! -f "$COLLAPSED" ]] || [[ "$COLLAPSED" -lt $(find results/*.vcf.gz -type f -printf '%T@ %p\n' | sort -n | tail -1 | cut -d' ' -f1) ]]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > "$COLLAPSED"
    for sample in "${SAMPLES[@]}"; do
        bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz | sed "s/^/${sample}\t/" >> "$COLLAPSED"
    done
fi