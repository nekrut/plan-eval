#!/usr/bin/env bash
set -euo pipefail

THREADS=4

# Ensure results directory exists
mkdir -p results

# Reference index
REF=data/ref/chrM.fa
if ! [[ -f ${REF}.bwt && -f ${REF}.pac && -f ${REF}.ann && -f ${REF}.amb && -f ${REF}.sa ]]; then
    bwa index -p ${REF} ${REF}
fi

# Sample list
samples=(M117-bl M117-ch M117C1-bl M117C1-ch)

for sample in "${samples[@]}"; do
    bam=results/${sample}.bam
    bai=results/${sample}.bam.bai
    vcf=results/${sample}.vcf.gz
    tbi=results/${sample}.vcf.gz.tbi

    # Mapping
    if [[ ! -f $bam ]]; then
        bwa mem -t $THREADS ${REF} data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
        samtools view -bS - | \
        samtools sort -o $bam
    fi

    # BAM index
    if [[ ! -f $bai ]]; then
        samtools index $bam
    fi

    # Variant calling
    if [[ ! -f $vcf ]]; then
        lofreq call -f ${REF} -o - $bam | \
        bcftools view -O z -o $vcf
    fi

    # VCF index
    if [[ ! -f $tbi ]]; then
        tabix -p vcf $vcf
    fi
done

# Collapsed table
collapsed=results/collapsed.tsv
if [[ ! -f $collapsed ]]; then
    printf "sample\tchrom\tpos\tref\talt\taf\n" > $collapsed
    for sample in "${samples[@]}"; do
        vcf=results/${sample}.vcf.gz
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' $vcf | \
        awk -v s=$sample '{print s"\t"$0}' >> $collapsed
    done
fi

exit 0