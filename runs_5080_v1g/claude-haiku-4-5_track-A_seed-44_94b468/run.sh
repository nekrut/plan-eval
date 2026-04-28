#!/usr/bin/env bash
set -euo pipefail

# Step 1: Set globals and prepare results directory
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
mkdir -p results

# Step 2: Reference indexing (once, in data/ref/)
if [[ ! -f data/ref/chrM.fa.fai ]]; then
    samtools faidx data/ref/chrM.fa
fi

if [[ ! -f data/ref/chrM.fa.amb ]]; then
    bwa index data/ref/chrM.fa
fi

# Step 3-5: Per-sample alignment, sorting, and indexing
for sample in "${SAMPLES[@]}"; do
    if [[ ! -f results/${sample}.bam.bai ]]; then
        # Alignment with bwa mem and pipe to samtools sort
        bwa mem -t ${THREADS} -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            data/ref/chrM.fa \
            data/raw/${sample}_1.fq.gz \
            data/raw/${sample}_2.fq.gz | \
        samtools sort -@ ${THREADS} -o results/${sample}.bam
        
        # Index BAM
        samtools index -@ ${THREADS} results/${sample}.bam
    fi
done

# Step 6-7: Variant calling with lofreq call-parallel and VCF compression/indexing
for sample in "${SAMPLES[@]}"; do
    if [[ ! -f results/${sample}.vcf.gz.tbi ]]; then
        # Call variants
        lofreq call-parallel --pp-threads ${THREADS} --verbose \
            --ref data/ref/chrM.fa --out results/${sample}.vcf \
            --sig --bonf \
            results/${sample}.bam
        
        # Compress VCF
        bgzip results/${sample}.vcf
        
        # Index compressed VCF
        tabix -p vcf results/${sample}.vcf.gz
    fi
done

# Step 8: Collapse step - aggregate variants from all samples
if [[ ! -f results/collapsed.tsv ]] || \
   [[ results/collapsed.tsv -ot results/M117-bl.vcf.gz ]] || \
   [[ results/collapsed.tsv -ot results/M117-ch.vcf.gz ]] || \
   [[ results/collapsed.tsv -ot results/M117C1-bl.vcf.gz ]] || \
   [[ results/collapsed.tsv -ot results/M117C1-ch.vcf.gz ]]; then
    
    # Create temporary file for collapse
    temp_collapse=$(mktemp)
    
    # Write header
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "${temp_collapse}"
    
    # Query each VCF and append
    for sample in "${SAMPLES[@]}"; do
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
            results/${sample}.vcf.gz >> "${temp_collapse}"
    done
    
    # Move to final location
    mv "${temp_collapse}" results/collapsed.tsv
fi

exit 0