#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF="data/ref/chrM.fa"
REF_IDX="data/ref/chrM.fa.bwt"
RESULTS_DIR="results"

mkdir -p "$RESULTS_DIR"

# 2. Reference indexing (once)
if [ ! -f "$REF_IDX" ]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# 3-7. Per-sample processing
for sample in "${SAMPLES[@]}"; do
    bam_file="$RESULTS_DIR/${sample}.bam"
    bai_file="$RESULTS_DIR/${sample}.bam.bai"
    vcf_file="$RESULTS_DIR/${sample}.vcf"
    vcf_gz_file="$RESULTS_DIR/${sample}.vcf.gz"
    vcf_tbi_file="$RESULTS_DIR/${sample}.vcf.gz.tbi"

    # Skip if BAM and BAI already exist
    if [ -f "$bai_file" ]; then
        continue
    fi

    # 3. Alignment
    # Note: \t in the -R argument must be literal backslash-t for bwa
    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" "$REF" \
        "data/raw/${sample}_1.fq.gz" "data/raw/${sample}_2.fq.gz" | \
    # 4. Sort BAM
    samtools sort -@ "$THREADS" -o "$bam_file" -

    # 5. Index BAM
    samtools index -@ "$THREADS" "$bam_file"

    # 6. Variant calling
    lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$vcf_file" "$bam_file"

    # 7. Compress and index VCF
    bgzip -c "$vcf_file" > "$vcf_gz_file"
    tabix -p vcf "$vcf_gz_file"
    rm -f "$vcf_file"
done

# 8. Collapse step
collapsed_file="$RESULTS_DIR/collapsed.tsv"
need_collapse=false

for sample in "${SAMPLES[@]}"; do
    vcf_gz_file="$RESULTS_DIR/${sample}.vcf.gz"
    if [ ! -f "$vcf_gz_file" ]; then
        need_collapse=true
        break
    fi
done

if [ ! -f "$collapsed_file" ]; then
    need_collapse=true
fi

if $need_collapse; then
    # Write header
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$collapsed_file"

    # Append data for each sample
    for sample in "${SAMPLES[@]}"; do
        vcf_gz_file="$RESULTS_DIR/${sample}.vcf.gz"
        if [ -f "$vcf_gz_file" ]; then
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$vcf_gz_file" >> "$collapsed_file"
        fi
    done
fi