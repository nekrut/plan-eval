#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF=data/ref/chrM.fa
RAW=data/raw
OUT=results

mkdir -p "$OUT"

if [ ! -f "${REF}.fai" ]; then
    samtools faidx "$REF"
fi

if [ ! -f "${REF}.bwt" ] || [ ! -f "${REF}.amb" ] || [ ! -f "${REF}.ann" ] || [ ! -f "${REF}.pac" ] || [ ! -f "${REF}.sa" ]; then
    bwa index "$REF"
fi

for sample in "${SAMPLES[@]}"; do
    bam="${OUT}/${sample}.bam"
    bai="${OUT}/${sample}.bam.bai"
    vcf="${OUT}/${sample}.vcf"
    vcfgz="${OUT}/${sample}.vcf.gz"
    tbi="${OUT}/${sample}.vcf.gz.tbi"
    fq1="${RAW}/${sample}_1.fq.gz"
    fq2="${RAW}/${sample}_2.fq.gz"

    if [ ! -f "$bam" ] || [ "$fq1" -nt "$bam" ] || [ "$fq2" -nt "$bam" ]; then
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            "$REF" "$fq1" "$fq2" \
            | samtools sort -@ "$THREADS" -o "$bam" -
    fi

    if [ ! -f "$bai" ] || [ "$bam" -nt "$bai" ]; then
        samtools index -@ "$THREADS" "$bam"
    fi

    if [ ! -f "$vcfgz" ] || [ "$bam" -nt "$vcfgz" ]; then
        rm -f "$vcf" "$vcfgz" "$tbi"
        lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$vcf" "$bam"
        bgzip "$vcf"
    fi

    if [ ! -f "$tbi" ] || [ "$vcfgz" -nt "$tbi" ]; then
        tabix -p vcf "$vcfgz"
    fi
done

collapsed="${OUT}/collapsed.tsv"
need_rebuild=0
if [ ! -f "$collapsed" ]; then
    need_rebuild=1
else
    for sample in "${SAMPLES[@]}"; do
        if [ "${OUT}/${sample}.vcf.gz" -nt "$collapsed" ]; then
            need_rebuild=1
            break
        fi
    done
fi

if [ "$need_rebuild" -eq 1 ]; then
    tmp="${collapsed}.tmp"
    {
        printf 'sample\tchrom\tpos\tref\talt\taf\n'
        for sample in "${SAMPLES[@]}"; do
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "${OUT}/${sample}.vcf.gz"
        done
    } > "$tmp"
    mv "$tmp" "$collapsed"
fi