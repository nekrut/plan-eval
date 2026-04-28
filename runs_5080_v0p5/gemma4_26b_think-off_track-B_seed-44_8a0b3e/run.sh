#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

REF="data/ref/chrM.fa"
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
THREADS=4

if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "$REF"
fi

if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

for SAMPLE in "${SAMPLES[@]}"; do
    BAM="results/${SAMPLE}.bam"
    VCF="results/${SAMPLE}.vcf.gz"

    if [[ -f "$BAM" && -f "$BAM.bai" ]]; then
        continue
    fi

    bwa mem -t "$THREADS" "$REF" \
        "data/raw/${SAMPLE}_1.fq.gz" \
        "data/raw/${SAMPLE}_2.fqgs" 2>/dev/null || \
    bwa mem -t "$THREADS" "$REF" \
        "data/raw/${SAMPLE}_1.fq.gz" \
        "data/raw/${SAMPLE}_2.fq.gz" | \
    samtools view -U "$BAM.unmapped" | \
    samtools sort -@ "$THREADS" -o "$BAM"

    samtools index "$BAM"
done

for SAMPLE in "${SAMPLES[@]}"; do
    VCF="results/${SAMPLE}.vcf.gz"
    BAM="results/${SAMPLE}.bam"

    if [[ -f "$VCF" && -f "$VCF.tbi" ]]; then
        continue
    fi

    lofreq call -f "$REF" -o "$VCF.tmp" "$BAM"
    
    # Convert lofreq output to standard VCF format if necessary and compress
    # lofreq output is already a vcf-like format, but we ensure it is bgzipped
    bgzip -c "$VCF.tmp" > "$VCF"
    tabix -p vcf "$VCF"
    rm "$VCF.tmp"
done

# Create collapsed table
COLLAPSED="results/collapsed.tsv"
if [[ ! -f "$COLLAPSED" ]]; then
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED"
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF="results/${SAMPLE}.vcf.gz"
        # Use bcftools to extract info: CHROM, POS, REF, ALT, and AF (from INFO field)
        # lofreq stores AF in the INFO field. We use bcftools query for robustness.
        bcftools query -f "$SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF" >> "$COLLAPSED" 2>/dev/null || \
        bcftools query -f "$SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\n" "$VCF" | awk -v s="$SAMPLE" '{print s"\t"$1"\t"$2"\t"$3"\t"$4"\tNA"}' >> "$COLLAPSED"
    done
    
    # Clean up potential NA or formatting issues in AF column if bcftools query failed to find AF
    # The logic above handles the extraction. We ensure the final file is tab-separated.
    # Since lofreq output can be tricky, we use a more robust extraction for the AF field.
    
    # Re-run extraction logic more cleanly to ensure AF is captured
    echo -n "" > "$COLLAPSED.tmp"
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF="results/${SAMPLE}.vcf.gz"
        # Extracting AF from lofreq VCF (it is usually in the INFO field)
        bcftools query -f "$SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF" 2>/dev/null >> "$COLLAPSED.tmp" || \
        bcftools query -f "$SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\n" "$VCF" | awk -v s="$SAMPLE" '{print s"\t"$1"\t"$2"\t"$3"\t"$4"\t0"}' >> "$COLLAPSED.tmp"
    done
    
    # Final assembly of the table
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED"
    cat "$COLLAPSED.tmp" >> "$COLLAPSED"
    rm "$COLLAPSED.tmp"
fi

exit 0