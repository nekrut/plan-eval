#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF="data/ref/chrM.fa"
RAW="data/raw"
OUT="results"
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

mkdir -p "$OUT"

if [ ! -f "${REF}.bwt" ]; then
    bwa index "$REF"
fi
if [ ! -f "${REF}.fai" ]; then
    samtools faidx "$REF"
fi

for sample in "${SAMPLES[@]}"; do
    bam="${OUT}/${sample}.bam"
    bai="${OUT}/${sample}.bam.bai"
    vcf="${OUT}/${sample}.vcf.gz"
    tbi="${OUT}/${sample}.vcf.gz.tbi"
    r1="${RAW}/${sample}_1.fq.gz"
    r2="${RAW}/${sample}_2.fq.gz"

    if [ ! -f "$bam" ] || [ ! -f "$bai" ]; then
        tmp_sorted="${OUT}/${sample}.sorted.tmp.bam"
        tmp_indel="${OUT}/${sample}.indel.tmp.bam"
        bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" "$REF" "$r1" "$r2" \
            | samtools sort -@ "$THREADS" -o "$tmp_sorted" -
        lofreq indelqual --dindel -f "$REF" -o "$tmp_indel" "$tmp_sorted"
        mv "$tmp_indel" "$bam"
        rm -f "$tmp_sorted"
        samtools index -@ "$THREADS" "$bam"
    fi

    if [ ! -f "$vcf" ] || [ ! -f "$tbi" ]; then
        tmp_vcf="${OUT}/${sample}.tmp.vcf"
        lofreq call-parallel --pp-threads "$THREADS" --call-indels -f "$REF" -o "$tmp_vcf" "$bam"
        bgzip -f "$tmp_vcf"
        mv "${tmp_vcf}.gz" "$vcf"
        tabix -p vcf "$vcf"
    fi
done

collapsed="${OUT}/collapsed.tsv"
if [ ! -f "$collapsed" ]; then
    tmp_collapsed="${collapsed}.tmp"
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$tmp_collapsed"
    for sample in "${SAMPLES[@]}"; do
        vcf="${OUT}/${sample}.vcf.gz"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$vcf" \
            | awk -v s="$sample" 'BEGIN{OFS="\t"} {print s,$1,$2,$3,$4,$5}' >> "$tmp_collapsed"
    done
    mv "$tmp_collapsed" "$collapsed"
fi