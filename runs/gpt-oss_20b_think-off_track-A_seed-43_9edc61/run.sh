#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# Reference indexing – BWA
[[ -f data/ref/chrM.fa.bwt ]] || bwa index data/ref/chrM.fa

# Reference indexing – samtools faidx
[[ -f data/ref/chrM.fa.fai ]] || samtools faidx data/ref/chrM.fa

for sample in "${SAMPLES[@]}"; do
    BAM="results/${sample}.bam"
    BAI="results/${sample}.bam.bai"
    VCF="results/${sample}.vcf"
    VCFZ="results/${sample}.vcf.gz"
    TBZ="results/${sample}.vcf.gz.tbi"

    # Alignment and sort
    if [[ ! -f "$BAM" ]]; then
        bwa mem -t "$THREADS" -R "@RG\\tID:${sample}\\tSM:${sample}\\tLB:${sample}\\tPL:ILLUMINA" \
            data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
        samtools sort -@ "$THREADS" -o "$BAM" -
    fi

    # BAM index
    [[ -f "$BAI" ]] || samtools index -@ "$THREADS" "$BAM"

    # Variant calling – LoFreq
    if [[ ! -f "$VCF" && ! -f "$VCFZ" ]]; then
        lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o "$VCF" "$BAM"
    fi

    # VCF compression and tabix index
    if [[ ! -f "$TBZ" ]]; then
        bgzip -f "$VCF"
        tabix -p vcf "$VCFZ"
    fi
done

# Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
for sample in "${SAMPLES[@]}"; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${sample}.vcf.gz" | \
    awk -v s="${sample}" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done