#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
mkdir -p results

# 2. Reference indexing (once)
if [[ ! -f data/ref/chrM.fa.fai ]]; then
    samtools faidx data/ref/chrM.fa
fi
if [[ ! -f data/ref/chrM.fa.bwt ]]; then
    bwa index data/ref/chrM.fa
fi

# 3. Per-sample alignment
for sample in "${SAMPLES[@]}"; do
    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
        data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
        samtools sort -@ "$THREADS" -o "results/${sample}.bam"
done

# 5. BAM indexing
for sample in "${SAMPLES[@]}"; do
    samtools index -@ "$THREADS" "results/${sample}.bam"
done

# 6. Variant calling
for sample in "${SAMPLES[@]}"; do
    lofreq call-parallel --pp-threads "$THREADS" --verbose \
        --ref data/ref/chrM.fa --out "results/${sample}.vcf" \
        --sig --bonf "results/${sample}.bam"
done

# 7. VCF compression and indexing
for sample in "${SAMPLES[@]}"; do
    bgzip -c "results/${sample}.vcf" > "results/${sample}.vcf.gz"
    tabix -p vcf "results/${sample}.vcf.gz"
    rm "results/${sample}.vcf"
done

# 8. Collapse step
if [[ ! -f results/collapsed.tsv ]] || \
   [[ results/collapsed.tsv -ot results/M117-bl.vcf.gz ]] || \
   [[ results/collapsed.tsv -ot results/M117-ch.vcf.gz ]] || \
   [[ results/collapsed.tsv -ot results/M117C1-bl.vcf.gz ]] || \
   [[ results/collapsed.tsv -ot results/M117C1-ch.vcf.gz ]]; then
    {
        echo -e "sample\tchrom\tpos\tref\talt\taf"
        for sample in "${SAMPLES[@]}"; do
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "results/${sample}.vcf.gz"
        done
    } > results/collapsed.tsv
fi