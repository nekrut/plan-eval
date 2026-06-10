#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

REF="data/ref/chrM.fa"
THREADS=4

# Index reference if not already indexed
if [ ! -f "${REF}._1.bam" ] && [ ! -f "${REF}.bt" ]; then
    # Note: bwa index requires .1.bam or similar, but standard is .amb, .ann, .pac, .sa
    # Since we use bwa mem and samtools, we need the standard indices.
    bwa index "$REF"
fi

# Create a dummy index for samtools if needed (though bwa index usually handles this)
if [ ! -f "${REF}.fai" ]; then
    samtools faidx "$REF"
fi

SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

for SAMPLE in "${SAMPLES[@]}"; do
    BAM="results/${SAMPLE}.bam"
    VCF="results/${SAMPLE}.vcf.gz"
    
    # Skip if already exists
    if [[ -f "$BAM" && -f "$VCF" ]]; then
        continue
    fi

    # Map and Sort/Index
    # Using bwa mem for alignment, samtools for processing
    bwa mem -t "$THREADS" -K "${REF}" "data/raw/${SAMPLE}_1.fq.gz" "data/raw/${SAMPLE}_2.fq.gz" | \
    samtools view -u - | \
    samtools sort -@ "$THREADS" -o "$BAM"
    
    samtools index "$BAM"

    # Variant Calling using lofreq (robust for low coverage/mitochondrial)
    # lofreq outputs a VCF. We then compress and index it.
    lofreq call -t "$THREADS" -f "$REF" "$BAM" > "${VCF}.tmp"
    mv "${VCF}.tmp" "$VCF"
    bgzip -c "$VCF" > "${VCF}.gz" # lofreq might not produce .gz, but we need it. 
    # Correction: lofreq output is usually a raw vcf. Let's ensure standard format.
    # Since the requirement asks for .vcf.gz and .tbi:
    
    # Re-running logic to ensure exact filenames match requirements
    # The previous step produced $VCF (which was actually a .tmp file)
    # Let's refine the flow to be strictly idempotent and correct.
done

# Refined loop for strict requirement adherence
for SAMPLE in "${SAMPLES[@]}"; do
    BAM="results/${SAMPLE}.bam"
    VCF_RAW="results/${SAMPLE}.vcf"
    VCF_FINAL="results/${SAMPLE}.vcf.gz"

    if [[ -f "$BAM" && -f "$VCF_FINAL" ]]; then
        continue
    fi

    # Alignment
    bwa mem -t "$THREADS" -K "$REF" "data/raw/${SAMPLE}_1.fq.gz" "data/raw/${SAMPLE}_2.fq.gz" | \
    samtools view -u - | \
    samtools sort -@ "$THREADS" -o "$BAM"
    samtools index "$BAM"

    # Calling
    lofreq call -t "$THREADS" -f "$REF" "$BAM" > "$VCF_RAW"
    
    # Post-processing VCF
    # lofreq output might have headers or specific formatting, but we need to compress it.
    # We use bcftools to ensure the VCF is standard before compression if needed.
    bcftools view -Oz "$VCF_RAW" > "$VCF_FINAL"
    rm "$VCF_RAW"
    tabix -1 "$VCF_FINAL"
done

# Create collapsed table
# Columns: sample, chrom, pos, ref, alt, af
# We use bcftools query to extract fields. 
# Note: lofreq provides AF in the INFO field or as a separate field depending on version.
# Standard VCF format for AF is usually in INFO.

echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv

for SAMPLE in "${SAMPLES[@]}"; do
    VCF="results/${SAMPLE}.vcf.gz"
    # Extract fields: SampleName, CHROM, POS, REF, ALT, and AF from INFO
    # If AF is not in INFO, we default to 1.0 or similar if available.
    bcftools query -mv -h "$VCF" | awk -v s="$SAMPLE" '{print s "\t" $1 "\t" $2 "\t" $3 "\t" $4}' >> results/collapsed.tsv.tmp
done

# Since bcftools query doesn't easily grab a specific INFO field like 'AF' 
# without complex formatting, we use a more robust approach to get the AF:
# We iterate through each VCF and extract the data.

# Clear the temp file and rebuild correctly
rm -f results/collapsed.tsv.tmp
echo -e "sample\tchrom\tpos\ref\talt\taf" > results/collapsed.tsv

for SAMPLE in "${SAMPLES[@]}"; do
    VCF="results/${SAMPLE}.vcf.gz"
    # Use bcftools query to get: Sample, Chrom, Pos, Ref, Alt, and the AF info field
    bcftools query -l "$VCF" | while read -r line; do :; done # dummy
    # Correct logic for extracting specific INFO fields:
    bcftools query -mv -p CHROM,POS,REF,ALT,AF "$VCF" | \
    awk -v s="$SAMPLE" '{print s "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' >> results/collapsed.tsv.tmp
done

# Final merge of the gathered data into the final file (excluding header from temp)
tail -n +2 results/collapsed.tsv.tmp >> results/collapsed.tsv
rm results/collapsed.tsv.tmp