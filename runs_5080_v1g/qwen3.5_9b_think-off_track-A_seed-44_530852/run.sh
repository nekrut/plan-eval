#!/usr/bin/env bash
set -euo pipefail
THREADS=4
SAMPLES="M117-bl M117-ch M117C1-bl M117C1-ch"
REF="data/ref/chrM.fa"
RAW="data/raw"
RES="results"

mkdir -p "$RES"

# 2. Reference indexing (once)
if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.amb" ]]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# 3-5. Per-sample alignment and BAM indexing
for sample in $SAMPLES; do
    f1="${RAW}/${sample}_1.fq.gz"
    f2="${RAW}/${sample}_2.fq.gz"
    out="${RES}/${sample}.bam"
    out_bai="${RES}/${sample}.bam.bai"
    
    if [[ ! -f "$out" ]]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" "$REF" "$f1" "$f2" | \
            samtools sort -@ "$THREADS" -o "$out"
        samtools index -@ "$THREADS" "$out"
    fi
done

# 6. Variant calling with lofreq
for sample in $SAMPLES; do
    out="${RES}/${sample}.vcf"
    if [[ ! -f "$out" ]]; then
        lofreq call-parallel --pp-threads "$THREADS" --verbose \
            --ref "$REF" --out "$out" "$RES/${sample}.bam"
    fi
done

# 7. VCF compression and indexing
for sample in $SAMPLES; do
    out="${RES}/${sample}.vcf"
    out_gz="${RES}/${sample}.vcf.gz"
    out_tbi="${RES}/${sample}.vcf.gz.tbi"
    
    if [[ ! -f "$out_gz" ]]; then
        bgzip -f "$out"
        tabix -p vcf "$out_gz"
        rm -f "$out"
    fi
done

# 8. Collapse step
TSV="${RES}/collapsed.tsv"
TMP="${RES}/collapsed.tmp"

{
    echo -e "sample\tchrom\tpos\tref\talt\taf"
    for sample in $SAMPLES; do
        bcftools query -f '%sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${RES}/${sample}.vcf.gz"
    done
} > "$TMP"

mv "$TMP" "$TSV"