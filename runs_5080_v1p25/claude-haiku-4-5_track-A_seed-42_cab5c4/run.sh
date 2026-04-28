#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF="data/ref/chrM.fa"
RESULTS_DIR="results"

mkdir -p "${RESULTS_DIR}"

# Reference indexing
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "${REF}"
fi

if [[ ! -f "${REF}.amb" ]]; then
    bwa index "${REF}"
fi

# Per-sample processing
for sample in "${SAMPLES[@]}"; do
    read1="data/raw/${sample}_1.fq.gz"
    read2="data/raw/${sample}_2.fq.gz"
    bam="${RESULTS_DIR}/${sample}.bam"
    bai="${bam}.bai"
    vcf="${RESULTS_DIR}/${sample}.vcf"
    vcf_gz="${RESULTS_DIR}/${sample}.vcf.gz"
    vcf_tbi="${vcf_gz}.tbi"
    
    # Alignment and sorting
    if [[ ! -f "${bam}" ]] || [[ "${read1}" -nt "${bam}" ]] || [[ "${read2}" -nt "${bam}" ]]; then
        bwa mem -t ${THREADS} -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" "${REF}" "${read1}" "${read2}" | \
            samtools sort -@ ${THREADS} -o "${bam}"
    fi
    
    # BAM indexing
    if [[ ! -f "${bai}" ]] || [[ "${bam}" -nt "${bai}" ]]; then
        samtools index -@ ${THREADS} "${bam}"
    fi
    
    # Variant calling, compression, and indexing
    if [[ ! -f "${vcf_gz}" ]] || [[ "${bam}" -nt "${vcf_gz}" ]]; then
        lofreq call-parallel --pp-threads ${THREADS} -f "${REF}" -o "${vcf}" "${bam}"
        rm -f "${vcf_gz}" "${vcf_tbi}"
        bgzip "${vcf}"
        tabix -p vcf "${vcf_gz}"
    fi
done

# Collapse step
collapsed="${RESULTS_DIR}/collapsed.tsv"

rebuild_collapsed=false
if [[ ! -f "${collapsed}" ]]; then
    rebuild_collapsed=true
else
    for sample in "${SAMPLES[@]}"; do
        vcf_gz="${RESULTS_DIR}/${sample}.vcf.gz"
        if [[ "${vcf_gz}" -nt "${collapsed}" ]]; then
            rebuild_collapsed=true
            break
        fi
    done
fi

if [[ "${rebuild_collapsed}" == true ]]; then
    {
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' sample chrom pos ref alt af
        for sample in "${SAMPLES[@]}"; do
            vcf_gz="${RESULTS_DIR}/${sample}.vcf.gz"
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "${vcf_gz}"
        done
    } > "${collapsed}"
fi

exit 0