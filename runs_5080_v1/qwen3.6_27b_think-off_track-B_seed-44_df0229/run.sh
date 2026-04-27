#!/usr/bin/env bash
set -euo pipefail

# Configuration
REF="data/ref/chrM.fa"
RAW_DIR="data/raw"
OUT_DIR="results"
THREADS=4

# Sample names (prefixes for R1/R2)
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Create output directory if it doesn't exist
mkdir -p "${OUT_DIR}"

# 1. Index Reference (if not already indexed)
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "${REF}"
fi
if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "${REF}"
fi

# 2. Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    R1="${RAW_DIR}/${SAMPLE}_1.fq.gz"
    R2="${RAW_DIR}/${SAMPLE}_2.fq.gz"
    BAM="${OUT_DIR}/${SAMPLE}.bam"
    BAI="${BAM}.bai"
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    TBI="${VCF_GZ}.tbi"

    # Skip if final VCF index exists (idempotency check)
    if [[ -f "${TBI}" ]]; then
        continue
    fi

    # Alignment
    if [[ ! -f "${BAM}" ]]; then
        bwa mem -t "${THREADS}" "${REF}" "${R1}" "${R2}" | \
            samtools view -b -@ "${THREADS}" -o "${BAM}"
    fi

    # Sorting and Indexing
    if [[ ! -f "${BAI}" ]]; then
        samtools sort -@ "${THREADS}" -o "${BAM}" "${BAM}"
        samtools index -@ "${THREADS}" "${BAM}"
    fi

    # Variant Calling with Lofreq
    # Lofreq outputs VCF to stdout or file. We pipe to bgzip to create .vcf.gz directly.
    if [[ ! -f "${VCF_GZ}" ]]; then
        lofreq call --call-indels -f "${REF}" -o - "${BAM}" | \
            bgzip -c > "${VCF_GZ}"
    fi

    # Index VCF
    if [[ ! -f "${TBI}" ]]; then
        tabix -p vcf "${VCF_GZ}"
    fi
done

# 3. Collapse VCFs into a single TSV
COLLAPSED="${OUT_DIR}/collapsed.tsv"

# Write header
printf "sample\tchrom\tpos\tref\talt\taf\n" > "${COLLAPSED}"

for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    
    # Extract INFO field for AF (Allele Frequency)
    # bcftools query format:
    # %CHROM %POS %REF %ALT %INFO/AF
    # Note: Lofreq AF is in INFO/AF. If multiple alts, AF is comma-separated.
    # We need to handle multi-ALT sites by splitting them if necessary, 
    # but typically for mtDNA amplicons, we might see simple SNPs/Indels.
    # To be safe and match the requested columns (singular alt/af), 
    # we will split multi-ALT records into separate rows.
    
    # Use bcftools query to get fields. 
    # If AF is missing, we might need to calculate or default, but Lofreq usually provides it.
    # Format: sample chrom pos ref alt af
    
    # We use a small awk script to handle multi-ALT splitting if needed.
    # bcftools query output:
    # CHROM POS REF ALT AF
    # Note: AF might be '.' if not available.
    
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${VCF_GZ}" | \
    awk -v sample="${SAMPLE}" '
    BEGIN { OFS="\t" }
    {
        chrom = $1
        pos = $2
        ref = $3
        alt_str = $4
        af_str = $5
        
        # Split ALT and AF by comma
        n_alt = split(alt_str, alts, ",")
        
        # Handle AF: if it is a single value, apply to all. If comma-separated, match index.
        if (af_str == ".") {
            af_val = "."
        } else {
            n_af = split(af_str, afs, ",")
        }
        
        for (i=1; i<=n_alt; i++) {
            if (af_str == ".") {
                current_af = "."
            } else if (n_af > 1) {
                current_af = afs[i]
            } else {
                current_af = afs[1]
            }
            print sample, chrom, pos, ref, alts[i], current_af
        }
    }' >> "${COLLAPSED}"
done