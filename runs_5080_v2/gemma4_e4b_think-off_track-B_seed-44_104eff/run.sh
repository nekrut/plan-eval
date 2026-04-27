#!/usr/bin/env bash
set -euo pipefail

SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF_FA="data/ref/chrM.fa"
RESULTS_DIR="results"
THREADS=4

mkdir -p "$RESULTS_DIR"

# 1. Index the reference genome (chrM.fa)
if [ ! -f "$RESULTS_DIR/chrM.fa.fai" ]; then
    samtools faidx "$REF_FA"
fi

# 2. Alignment, Sorting, and Indexing for each sample
for sample in "${SAMPLES[@]}"; do
    echo "Processing sample: $sample"

    # Check if BAM files already exist (idempotency check)
    if [ -f "$RESULTS_DIR/${sample}.bam" ]; then
        echo "Skipping alignment/sorting for $sample: BAM already exists."
        continue
    fi

    # Concatenate paired-end reads
    R1="${sample}_1.fq.gz"
    R2="${sample}_2.fq.gz"

    if [ ! -f "data/raw/${R1}" ] || [ ! -f "data/raw/${R2}" ]; then
        echo "Error: Missing input files for $sample. Skipping."
        continue
    fi

    # BWA MEM mapping
    echo "Mapping reads for $sample..."
    bwa mem -t $THREADS -R "@RG\tID:$sample\tSM:$sample\tPL:ILLUMINA" "$REF_FA" "data/raw/${R1}" "data/raw/${R2}" | \
    samtools view -@ $THREADS -Sb - | \
    samtools sort -@ $THREADS -o "$RESULTS_DIR/${sample}.sorted.bam" -

    # Indexing
    samtools index "$RESULTS_DIR/${sample}.sorted.bam"
    mv "$RESULTS_DIR/${sample}.sorted.bam" "$RESULTS_DIR/${sample}.bam"

    # Cleanup intermediate sorted file if necessary (though mv handles it)
    rm -f "$RESULTS_DIR/${sample}.sorted.bam"

    # Create VCF/BCF files
    if [ ! -f "$RESULTS_DIR/${sample}.vcf.gz" ]; then
        echo "Calling variants for $sample..."
        bcftools mpileup -f "$REF_FA" -g -t "$sample" -D "$RESULTS_DIR/${sample}.bam" | \
        bcftools call -mv -Oz -o "$RESULTS_DIR/${sample}.vcf.gz"
        tabix -p vcf "$RESULTS_DIR/${sample}.vcf.gz"
    fi
done

# 3. Collapse VCFs into a single TSV table
echo "Collapsing results into collapsed.tsv..."

# Initialize the collapsed file with header
echo -e "sample\tchrom\tpos\tref\talt\taf" > "$RESULTS_DIR/collapsed.tsv"

# Process each sample's VCF
for sample in "${SAMPLES[@]}"; do
    VCF="$RESULTS_DIR/${sample}.vcf.gz"
    if [ -f "$VCF" ]; then
        # Extract required fields: sample name, CHROM, POS, REF, ALT, AF
        # bcftools query is used to extract fields, then awk processes them into TSV format
        bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%AD\n" "$VCF" | \
        awk -v sample="$sample" '{
            # $1=CHROM, $2=POS, $3=REF, $4=ALT, $5=AD (read depth proxy for AF)
            # We need to handle multiple ALT alleles and calculate a representative AF.
            # For simplicity and adherence to the required columns, we will use the first ALT allele
            # and the sum of read depths (AD) as a proxy for 'af' if multiple are present.
            # Since the required output is a single TSV, we process one record at a time.
            
            # Reconstructing the required format: sample | chrom | pos | ref | alt | af
            # $1=CHROM, $2=POS, $3=REF, $4=ALT (comma-separated), $5=AD (read depth)
            
            # For this specific requirement, we assume the first ALT allele and use the total depth.
            # If ALT is complex, this will only capture the first part.
            alt_alleles = $4;
            if (alt_alleles ~ /,/) {
                split(alt_alleles, alts, ",");
                alt_out = alts[1]; # Take the first ALT allele
            } else {
                alt_out = alt_alleles;
            }
            
            # Use the sum of AD as a simple 'af' proxy if available, otherwise 0.
            af_out = $5;
            if (af_out == "") {
                af_out = 0;
            }
            
            print sample "\t" $1 "\t" $2 "\t" $3 "\t" alt_out "\t" af_out;
        }' >> "$RESULTS_DIR/collapsed.tsv"
    fi
done

# Final cleanup/verification (optional, but good practice)
# Ensure the header is correct and the file exists
if [ ! -f "$RESULTS_DIR/collapsed.tsv" ]; then
    echo "Warning: collapsed.tsv was not created."
fi

exit 0