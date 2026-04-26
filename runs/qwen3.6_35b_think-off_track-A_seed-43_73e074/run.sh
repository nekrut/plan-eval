#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF="data/ref/chrM.fa"
REF_DIR="data/ref"
RAW_DIR="data/raw"
RES_DIR="results"

mkdir -p "$RES_DIR"

# 2. Reference indexing (once)
if [ ! -f "$REF_DIR/chrM.fa.fai" ]; then
    samtools faidx "$REF"
fi
if [ ! -f "$REF_DIR/chrM.fa.bwt" ]; then
    bwa index "$REF"
fi

# 3-7. Per-sample processing
for sample in "${SAMPLES[@]}"; do
    bam="$RES_DIR/${sample}.bam"
    bai="$RES_DIR/${sample}.bam.bai"
    vcf_gz="$RES_DIR/${sample}.vcf.gz"
    vcf_tbi="$RES_DIR/${sample}.vcf.gz.tbi"
    vcf="$RES_DIR/${sample}.vcf"

    # 3-4. Alignment and sorting
    if [ ! -f "$bam" ] || [ "$RAW_DIR/${sample}_1.fq.gz" -nt "$bam" ] || [ "$RAW_DIR/${sample}_2.fq.gz" -nt "$bam" ]; then
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            "$REF" \
            "$RAW_DIR/${sample}_1.fq.gz" \
            "$RAW_DIR/${sample}_2.fq.gz" \
            | samtools sort -@ "$THREADS" -o "$bam" -
    fi

    # 5. BAM indexing
    if [ ! -f "$bai" ] || [ "$bam" -nt "$bai" ]; then
        samtools index -@ "$THREADS" "$bam"
    fi

    # 6. Variant calling with lofreq
    if [ ! -f "$vcf" ] || [ "$bam" -nt "$vcf" ]; then
        lofreq call-parallel -f "$REF" -i "$bam" -o "$vcf" --pp-threads "$THREADS"
    fi

    # 7. VCF compression and indexing
    if [ ! -f "$vcf_gz" ] || [ "$vcf" -nt "$vcf_gz" ]; then
        bgzip -c "$vcf" > "$vcf_gz"
        rm -f "$vcf"
    fi
    if [ ! -f "$vcf_tbi" ] || [ "$vcf_gz" -nt "$vcf_tbi" ]; then
        tabix -p vcf "$vcf_gz"
    fi
done

# 8. Collapse step
collapsed="$RES_DIR/collapsed.tsv"
need_collapse=false
for sample in "${SAMPLES[@]}"; do
    if [ ! -f "$collapsed" ] || [ "$RES_DIR/${sample}.vcf.gz" -nt "$collapsed" ]; then
        need_collapse=true
        break
    fi
done

if [ "$need_collapse" = true ]; then
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$collapsed"
    for sample in "${SAMPLES[@]}"; do
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$RES_DIR/${sample}.vcf.gz" >> "$collapsed"
    done
fi