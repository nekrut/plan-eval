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

    if [ ! -s "$bam" ] || [ ! -s "$bai" ]; then
        tmp_bam="${OUT}/${sample}.unviterbi.bam"
        rg="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
        bwa mem -t "$THREADS" -R "$rg" "$REF" "$r1" "$r2" \
            | samtools sort -@ "$THREADS" -o "$tmp_bam" -
        viterbi_bam="${OUT}/${sample}.indelqual.bam"
        lofreq viterbi -f "$REF" "$tmp_bam" \
            | samtools sort -@ "$THREADS" -o "${OUT}/${sample}.viterbi.bam" -
        lofreq indelqual --dindel -f "$REF" -o "$viterbi_bam" "${OUT}/${sample}.viterbi.bam"
        mv "$viterbi_bam" "$bam"
        samtools index -@ "$THREADS" "$bam"
        rm -f "$tmp_bam" "${OUT}/${sample}.viterbi.bam"
    fi

    if [ ! -s "$vcf" ] || [ ! -s "$tbi" ]; then
        raw_vcf="${OUT}/${sample}.raw.vcf"
        lofreq call-parallel --pp-threads "$THREADS" \
            --call-indels \
            -f "$REF" \
            -o "$raw_vcf" \
            "$bam"
        bcftools view -Oz -o "$vcf" "$raw_vcf"
        tabix -p vcf "$vcf"
        rm -f "$raw_vcf"
    fi
done

collapsed="${OUT}/collapsed.tsv"
if [ ! -s "$collapsed" ]; then
    tmp_collapsed="${collapsed}.tmp"
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$tmp_collapsed"
    for sample in "${SAMPLES[@]}"; do
        vcf="${OUT}/${sample}.vcf.gz"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$vcf" \
            | awk -v s="$sample" 'BEGIN{OFS="\t"} {print s,$1,$2,$3,$4,$5}' \
            >> "$tmp_collapsed"
    done
    mv "$tmp_collapsed" "$collapsed"
fi