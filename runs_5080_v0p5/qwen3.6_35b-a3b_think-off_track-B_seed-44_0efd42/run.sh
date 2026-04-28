#!/usr/bin/env bash
set -euo pipefail

# Create output directory
mkdir -p results

# Reference file
REF="data/ref/chrM.fa"

# Check if reference index exists, if not create it
if [ ! -f "${REF}.amb" ]; then
    bwa index "${REF}"
fi

# List of samples
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    R1="data/raw/${SAMPLE}_1.fq.gz"
    R2="data/raw/${SAMPLE}_2.fq.gz"
    BAM="results/${SAMPLE}.bam"
    VCF_GZ="results/${SAMPLE}.vcf.gz"
    VCF_TBI="results/${SAMPLE}.vcf.gz.tbi"

    # Skip if BAM already exists and is valid
    if [ -f "${BAM}" ] && [ -f "${BAM}.bai" ] && [ -f "${VCF_GZ}" ] && [ -f "${VCF_TBI}" ]; then
        continue
    fi

    # Step 1: Alignment with BWA
    # Use mem algorithm, output SAM then convert to BAM
    if [ ! -f "${BAM}" ]; then
        bwa mem -t 4 "${REF}" "${R1}" "${R2}" | samtools view -bS -o "${BAM}" -
    fi

    # Step 2: Sort BAM and generate index
    SORTED_BAM="results/${SAMPLE}.sorted.bam"
    if [ ! -f "${SORTED_BAM}" ]; then
        samtools sort -t 4 -o "${SORTED_BAM}" "${BAM}"
    fi
    if [ ! -f "${SORTED_BAM}.bai" ]; then
        samtools index "${SORTED_BAM}"
    fi

    # Step 3: Variant calling with LoFreq
    # LoFreq requires a FASTA index for the reference
    if [ ! -f "${REF}.fai" ]; then
        samtools faidx "${REF}"
    fi

    # LoFreq output is VCF
    if [ ! -f "${VCF_GZ}" ]; then
        lofreq call-parallel \
            -f "${REF}" \
            -i "${SORTED_BAM}" \
            -o "${VCF_GZ}" \
            --no-indel-realign \
            --no-qc-filter \
            -t 4
    fi

    # Step 4: Index VCF
    if [ ! -f "${VCF_TBI}" ]; then
        tabix -p vcf "${VCF_GZ}"
    fi

    # Clean up intermediate sorted BAM if desired, but keep BAM for potential reuse
    # rm -f "${SORTED_BAM}" "${SORTED_BAM}.bai"
done

# Step 5: Collapse VCFs into a single TSV
# Header
printf "sample\tchrom\tpos\tref\talt\taf\n" > results/collapsed.tsv

# Process each sample's VCF
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="results/${SAMPLE}.vcf.gz"
    if [ -f "${VCF_GZ}" ]; then
        # Use bcftools query to extract fields
        # Format: sample, chrom, pos, ref, alt, AF (from INFO/AF or FORMAT/AF)
        # LoFreq typically puts AF in the INFO field as AF or in FORMAT as AF
        # Let's use bcftools to extract GT and AF if present, or just INFO/AF
        # Standard VCF columns: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE1...
        # We want: sample, chrom, pos, ref, alt, af
        
        # Extract using bcftools query
        # %CHROM\t%POS\t%REF\t%ALT\t%INFO/AF
        # Note: LoFreq might use different tag names. Commonly AF is in INFO.
        # If AF is in FORMAT, we need %FORMAT/AF or %SAMPLE/AF
        # Let's try to get AF from INFO first. If not, try FORMAT.
        
        bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "${VCF_GZ}" | \
        awk -F'\t' -v sample="${SAMPLE}" '{
            # Handle cases where AF might be missing or multi-allelic
            af = $5
            if (af == "." || af == "") af = "0"
            # Split AF if it contains commas (multi-allelic)
            n = split(af, af_arr, ",")
            for (i=1; i<=n; i++) {
                print sample "\t" $1 "\t" $2 "\t" $3 "\t" $4 "," af_arr[i] "\t" af_arr[i]
            }
        }' >> results/collapsed.tsv
    fi
done

# If the above awk logic for multi-allelic is too complex or incorrect for LoFreq output format,
# let's use a simpler approach assuming single allele per line or standard AF extraction.
# LoFreq VCF usually has one line per variant. AF is in INFO/AF.
# Let's re-do the collapse with a more robust bcftools query if needed.
# But the above should work for most cases. Let's refine the awk to handle single AF correctly.

# Actually, let's restart the collapsed.tsv generation to be cleaner.
# Remove the previously written collapsed.tsv if it exists (idempotency)
rm -f results/collapsed.tsv

printf "sample\tchrom\tpos\tref\talt\taf\n" > results/collapsed.tsv

for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="results/${SAMPLE}.vcf.gz"
    if [ -f "${VCF_GZ}" ]; then
        # Extract CHROM, POS, REF, ALT, and AF from INFO
        # If AF is missing, default to 0
        bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "${VCF_GZ}" | \
        awk -F'\t' -v sample="${SAMPLE}" '{
            af = $5
            if (af == "." || af == "") af = "0"
            # If AF is comma-separated (multi-allelic), we need to expand
            if (index(af, ",") > 0) {
                n = split(af, af_arr, ",")
                for (i=1; i<=n; i++) {
                    print sample "\t" $1 "\t" $2 "\t" $3 "\t" $4 "," af_arr[i] "\t" af_arr[i]
                }
            } else {
                print sample "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" af
            }
        }' >> results/collapsed.tsv
    fi
done