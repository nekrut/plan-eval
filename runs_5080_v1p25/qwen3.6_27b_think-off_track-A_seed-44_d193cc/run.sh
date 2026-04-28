#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF="data/ref/chrM.fa"
RAW_DIR="data/raw"
OUT_DIR="results"

mkdir -p "$OUT_DIR"

# 2. Reference indexing (idempotent)
if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.bwt" ]]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# 3-7. Per-sample alignment, sorting, indexing, calling, compression
for sample in "${SAMPLES[@]}"; do
    bam_file="${OUT_DIR}/${sample}.bam"
    bai_file="${OUT_DIR}/${sample}.bam.bai"
    vcf_gz_file="${OUT_DIR}/${sample}.vcf.gz"
    tbi_file="${OUT_DIR}/${sample}.vcf.gz.tbi"

    # Idempotency guard: if final VCF index exists, skip sample
    if [[ -f "$tbi_file" ]]; then
        continue
    fi

    # 3 & 4. Alignment and sorting
    if [[ ! -f "$bam_file" ]]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            "$REF" \
            "${RAW_DIR}/${sample}_1.fq.gz" \
            "${RAW_DIR}/${sample}_2.fq.gz" | \
        samtools sort -@ "$THREADS" -o "$bam_file" -
    fi

    # 5. BAM indexing
    if [[ ! -f "$bai_file" ]]; then
        samtools index -@ "$THREADS" "$bam_file"
    fi

    # 6. Variant calling
    vcf_file="${OUT_DIR}/${sample}.vcf"
    if [[ ! -f "$vcf_file" ]]; then
        lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$vcf_file" "$bam_file"
    fi

    # 7. VCF compression and indexing
    if [[ ! -f "$vcf_gz_file" ]]; then
        bgzip -c "$vcf_file" > "$vcf_gz_file"
        tabix -p vcf "$vcf_gz_file"
        rm -f "$vcf_file"
    fi
done

# 8. Collapse step
COLLAPSED_TSV="${OUT_DIR}/collapsed.tsv"

# Check if any VCF is newer than the collapsed TSV
needs_collapse=false
if [[ ! -f "$COLLAPSED_TSV" ]]; then
    needs_collapse=true
else
    for sample in "${SAMPLES[@]}"; do
        vcf_gz_file="${OUT_DIR}/${sample}.vcf.gz"
        if [[ "$vcf_gz_file" -nt "$COLLAPSED_TSV" ]]; then
            needs_collapse=true
            break
        fi
    done
fi

if [[ "$needs_collapse" == true ]]; then
    {
        printf "sample\tchrom\tpos\tref\talt\taf\n"
        for sample in "${SAMPLES[@]}"; do
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "${OUT_DIR}/${sample}.vcf.gz"
        done
    } > "$COLLAPSED_TSV"
fi