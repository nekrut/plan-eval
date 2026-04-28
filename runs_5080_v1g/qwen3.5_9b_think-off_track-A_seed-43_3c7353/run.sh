#!/usr/bin/env bash
set -euo pipefail
THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF=data/ref/chrM.fa
RAW=data/raw
RES=results

mkdir -p "$RES"

# Reference indexing (idempotent)
if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.amb" ]]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# Per-sample processing
for SAMPLE in "${SAMPLES[@]}"; do
    F1="${RAW}/${SAMPLE}_1.fq.gz"
    F2="${RAW}/${SAMPLE}_2.fq.gz"
    OUT="${RES}/${SAMPLE}.bam"
    OUT_BAI="${OUT}.bai"
    OUT_VCF="${RES}/${SAMPLE}.vcf.gz"
    OUT_VCF_TBI="${OUT_VCF}.tbi"

    # Alignment (idempotent: skip if BAM exists)
    if [[ ! -f "$OUT" ]]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" "$REF" "$F1" "$F2" | \
            samtools sort -@ "$THREADS" -o "$OUT"
        samtools index -@ "$THREADS" "$OUT"
    fi

    # Variant calling (idempotent: skip if VCF exists)
    if [[ ! -f "$OUT_VCF" ]]; then
        lofreq call-parallel --pp-threads "$THREADS" --verbose \
            --ref "$REF" --out "$RES/${SAMPLE}.vcf" "$OUT"
        rm -f "${RES}/${SAMPLE}.vcf"
        bgzip "$RES/${SAMPLE}.vcf"
        tabix -p vcf "$OUT_VCF"
    fi
done

# Collapse step (idempotent: rebuild only if any VCF is newer)
TSV="${RES}/collapsed.tsv"
HEADER="sample	chrom	pos	ref	alt	af"
if [[ ! -f "$TSV" ]] || \
   [[ $(find "$RES" -name "*.vcf.gz" -newer "$TSV" | wc -l) -gt 0 ]]; then
    {
        echo -e "$HEADER"
        for SAMPLE in "${SAMPLES[@]}"; do
            bcftools query -f '%sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${RES}/${SAMPLE}.vcf.gz"
        done
    } > "$TSV"
fi