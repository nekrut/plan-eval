#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF=data/ref/chrM.fa
RAWDIR=data/raw
OUTDIR=results

SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

mkdir -p "${OUTDIR}"

# Index reference genome if not already done
if [ ! -f "${REF}.bwt" ]; then
    bwa index "${REF}"
fi

if [ ! -f "${REF}.fai" ]; then
    samtools faidx "${REF}"
fi

for SAMPLE in "${SAMPLES[@]}"; do
    BAM="${OUTDIR}/${SAMPLE}.bam"
    IQBAM="${OUTDIR}/${SAMPLE}.indelqual.bam"
    VCFGZ="${OUTDIR}/${SAMPLE}.vcf.gz"

    # Align, sort, and produce final BAM
    if [ ! -f "${BAM}" ]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
            "${REF}" \
            "${RAWDIR}/${SAMPLE}_1.fq.gz" \
            "${RAWDIR}/${SAMPLE}_2.fq.gz" \
        | samtools sort -@ "${THREADS}" -o "${BAM}" -
    fi

    # Index the final BAM
    if [ ! -f "${BAM}.bai" ]; then
        samtools index "${BAM}"
    fi

    # Add indel qualities required by lofreq
    if [ ! -f "${IQBAM}" ]; then
        lofreq indelqual --dindel -f "${REF}" -o "${IQBAM}" "${BAM}"
    fi

    if [ ! -f "${IQBAM}.bai" ]; then
        samtools index "${IQBAM}"
    fi

    # Variant calling with lofreq
    if [ ! -f "${VCFGZ}" ]; then
        RAWVCF="${OUTDIR}/${SAMPLE}.raw.vcf"
        lofreq call \
            --call-indels \
            -f "${REF}" \
            -o "${RAWVCF}" \
            "${IQBAM}"
        bgzip -c "${RAWVCF}" > "${VCFGZ}"
        rm -f "${RAWVCF}"
    fi

    # Index the compressed VCF
    if [ ! -f "${VCFGZ}.tbi" ]; then
        tabix -p vcf "${VCFGZ}"
    fi
done

# Build collapsed table across all samples
COLLAPSED="${OUTDIR}/collapsed.tsv"
if [ ! -f "${COLLAPSED}" ]; then
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${COLLAPSED}"
    for SAMPLE in "${SAMPLES[@]}"; do
        bcftools query \
            -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
            "${OUTDIR}/${SAMPLE}.vcf.gz" >> "${COLLAPSED}"
    done
fi