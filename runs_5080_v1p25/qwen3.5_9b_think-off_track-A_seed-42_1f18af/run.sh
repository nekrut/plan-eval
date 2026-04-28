#!/usr/bin/env bash
set -euo pipefail
THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF=data/ref/chrM.fa
RESULTS=results

mkdir -p "$RESULTS"

# Reference indexing (idempotent)
if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.amb" ]]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# Per-sample processing
for sample in "${SAMPLES[@]}"; do
    R1="data/raw/${sample}_1.fq.gz"
    R2="data/raw/${sample}_2.fq.gz"
    BAM="${RESULTS}/${sample}.bam"
    BAI="${RESULTS}/${sample}.bam.bai"
    VCF="${RESULTS}/${sample}.vcf.gz"
    VCF_TBI="${RESULTS}/${sample}.vcf.gz.tbi"
    
    # Alignment (idempotent: skip if BAM exists)
    if [[ ! -f "$BAM" ]]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" "$REF" "$R1" "$R2" 2>/dev/null | \
            samtools sort -@ "$THREADS" -o "$BAM"
        samtools index -@ "$THREADS" "$BAM"
    fi
    
    # Variant calling (idempotent: skip if VCF.gz.tbi exists)
    if [[ ! -f "$VCF_TBI" ]]; then
        lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "${RESULTS}/${sample}.vcf" "$BAM"
        bgzip -f "${RESULTS}/${sample}.vcf"
        tabix -p vcf "${RESULTS}/${sample}.vcf.gz"
        rm -f "${RESULTS}/${sample}.vcf"
    fi
done

# Collapse step (idempotent: rebuild only if any VCF is newer)
TSV="${RESULTS}/collapsed.tsv"
HEADER="sample	chrom	pos	ref	alt	af"
if [[ ! -f "$TSV" ]] || \
   [[ $(find "$RESULTS" -name "*.vcf.gz" -newer "$TSV" | wc -l) -gt 0 ]]; then
    {
        echo -e "$HEADER"
        for sample in "${SAMPLES[@]}"; do
            bcftools query -f '%sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${RESULTS}/${sample}.vcf.gz"
        done
    } > "$TSV"
fi