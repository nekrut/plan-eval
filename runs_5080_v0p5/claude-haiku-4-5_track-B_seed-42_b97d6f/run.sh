#!/usr/bin/env bash
set -euo pipefail

THREADS=${THREADS:-4}
mkdir -p results

REF="data/ref/chrM.fa"

# Index reference
if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "${REF}"
fi

SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    R1="data/raw/${SAMPLE}_1.fq.gz"
    R2="data/raw/${SAMPLE}_2.fq.gz"
    
    BAM="results/${SAMPLE}.bam"
    BAI="results/${SAMPLE}.bam.bai"
    VCF_TMP="results/${SAMPLE}.vcf"
    VCF="results/${SAMPLE}.vcf.gz"
    TBI="results/${SAMPLE}.vcf.gz.tbi"
    
    # Alignment and BAM creation
    if [[ ! -f "${BAM}" ]]; then
        bwa mem -t "${THREADS}" "${REF}" "${R1}" "${R2}" | \
            samtools view -b -h - | \
            samtools sort -o "${BAM}" -
    fi
    
    # Index BAM
    if [[ ! -f "${BAI}" ]]; then
        samtools index "${BAM}"
    fi
    
    # Variant calling
    if [[ ! -f "${VCF}" ]]; then
        lofreq call -f "${REF}" -o "${VCF_TMP}" "${BAM}"
        bgzip -f "${VCF_TMP}"
    fi
    
    # Index VCF
    if [[ ! -f "${TBI}" ]]; then
        tabix -p vcf "${VCF}"
    fi
done

# Collapse VCFs into single TSV
COLLAPSED="results/collapsed.tsv"
if [[ ! -f "${COLLAPSED}" ]]; then
    {
        echo -e "sample\tchrom\tpos\tref\talt\taf"
        for SAMPLE in "${SAMPLES[@]}"; do
            VCF="results/${SAMPLE}.vcf.gz"
            zcat "${VCF}" | awk -v s="${SAMPLE}" -F'\t' '!/^#/ {
                af = 0
                split($8, info_fields, ";")
                for (i = 1; i <= length(info_fields); i++) {
                    field = info_fields[i]
                    if (substr(field, 1, 3) == "AF=") {
                        af = substr(field, 4)
                        break
                    }
                }
                print s "\t" $1 "\t" $2 "\t" $4 "\t" $5 "\t" af
            }'
        done
    } > "${COLLAPSED}"
fi