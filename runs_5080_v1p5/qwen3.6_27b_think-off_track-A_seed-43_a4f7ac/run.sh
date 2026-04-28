#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
mkdir -p results

# 1. Reference indexing — BWA (idempotent)
if [[ ! -f data/ref/chrM.fa.bwt ]]; then
    bwa index data/ref/chrM.fa
fi

# 2. Reference indexing — samtools faidx (idempotent)
if [[ ! -f data/ref/chrM.fa.fai ]]; then
    samtools faidx data/ref/chrM.fa
fi

# 3-6. Per-sample alignment, sorting, indexing, calling, compression
for sample in "${SAMPLES[@]}"; do
    bam_out="results/${sample}.bam"
    bai_out="results/${sample}.bam.bai"
    vcf_out="results/${sample}.vcf.gz"
    tbi_out="results/${sample}.vcf.gz.tbi"

    # 3. Alignment + sort (idempotent: skip if BAM exists)
    if [[ ! -f "${bam_out}" ]]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            data/ref/chrM.fa \
            "data/raw/${sample}_1.fq.gz" \
            "data/raw/${sample}_2.fq.gz" \
        | samtools sort -@ "${THREADS}" -o "${bam_out}" -
    fi

    # 4. BAM index (idempotent: skip if BAI exists)
    if [[ ! -f "${bai_out}" ]]; then
        samtools index -@ "${THREADS}" "${bam_out}"
    fi

    # 5. Variant calling — LoFreq (idempotent: skip if VCF.gz exists)
    # Note: lofreq call-parallel outputs to stdout or file. We use -o to write VCF.
    # We check for the final compressed VCF to determine if calling is done.
    if [[ ! -f "${vcf_out}" ]]; then
        # lofreq call-parallel writes VCF to -o file.
        # We need to ensure the intermediate .vcf is handled.
        # lofreq call-parallel --pp-threads 4 -f ref -o out.vcf in.bam
        lofreq call-parallel --pp-threads "${THREADS}" -f data/ref/chrM.fa -o "results/${sample}.vcf" "${bam_out}"
        
        # 6. VCF compression + tabix index
        bgzip -f "results/${sample}.vcf"
        tabix -p vcf "${vcf_out}"
    fi
done

# 7. Collapsed TSV (idempotent: skip if exists)
if [[ ! -f results/collapsed.tsv ]]; then
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
    for sample in "${SAMPLES[@]}"; do
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${sample}.vcf.gz" | \
        awk -v s="${sample}" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
    done
fi