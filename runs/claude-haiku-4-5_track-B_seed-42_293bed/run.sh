#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF="data/ref/chrM.fa"
RESULTS="results"
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p "$RESULTS"

# Index reference
if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "$REF"
fi
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

# Process each sample
for sample in "${SAMPLES[@]}"; do
    BAM="$RESULTS/${sample}.bam"
    VCFGZ="$RESULTS/${sample}.vcf.gz"
    
    # Alignment
    if [[ ! -f "$BAM" ]]; then
        bwa mem -t "$THREADS" "$REF" \
            "data/raw/${sample}_1.fq.gz" \
            "data/raw/${sample}_2.fq.gz" | \
        samtools sort -@ "$THREADS" -o "$BAM"
        samtools index "$BAM"
    fi
    
    # Variant calling
    if [[ ! -f "$VCFGZ" ]]; then
        bcftools mpileup -f "$REF" "$BAM" | \
        bcftools call -m -v | \
        bcftools view -O z -o "$VCFGZ"
        tabix -p vcf "$VCFGZ"
    fi
done

# Collapse variants
if [[ ! -f "$RESULTS/collapsed.tsv" ]]; then
    {
        echo -e "sample\tchrom\tpos\tref\talt\taf"
        for sample in "${SAMPLES[@]}"; do
            VCFGZ="$RESULTS/${sample}.vcf.gz"
            if [[ -f "$VCFGZ" ]]; then
                bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n' "$VCFGZ" | \
                awk -v s="$sample" '{
                    n_alt = split($4, alt, ",")
                    n_ac = split($5, ac, ",")
                    an=$6
                    for (i=1; i<=n_alt; i++) {
                        af=0
                        if(an>0) af=ac[i]/an
                        print s"\t"$1"\t"$2"\t"$3"\t"alt[i]"\t"af
                    }
                }'
            fi
        done
    } > "$RESULTS/collapsed.tsv"
fi