#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF="data/ref/chrM.fa"
RESULTS="results"

mkdir -p "$RESULTS"

# Index reference
[[ -f "${REF}.bwt" ]] || bwa index "$REF"
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

# Process each sample
for SAMPLE in M117-bl M117-ch M117C1-bl M117C1-ch; do
    BAM="${RESULTS}/${SAMPLE}.bam"
    VCF="${RESULTS}/${SAMPLE}.vcf"
    VCF_GZ="${VCF}.gz"
    
    # Alignment
    if [[ ! -f "$BAM" ]]; then
        bwa mem -t "$THREADS" "$REF" \
            "data/raw/${SAMPLE}_1.fq.gz" \
            "data/raw/${SAMPLE}_2.fq.gz" | \
        samtools sort -@ "$THREADS" -o "$BAM" -
    fi
    
    # Index BAM
    [[ -f "${BAM}.bai" ]] || samtools index "$BAM"
    
    # Variant calling
    if [[ ! -f "$VCF_GZ" ]]; then
        [[ -f "$VCF" ]] && rm "$VCF"
        lofreq call -f "$REF" -o "$VCF" "$BAM"
        bgzip "$VCF"
    fi
    
    # Index VCF
    [[ -f "${VCF_GZ}.tbi" ]] || tabix -p vcf "$VCF_GZ"
done

# Generate collapsed.tsv
COLLAPSED="${RESULTS}/collapsed.tsv"
if [[ ! -f "$COLLAPSED" ]]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > "$COLLAPSED"
    for SAMPLE in M117-bl M117-ch M117C1-bl M117C1-ch; do
        zcat "${RESULTS}/${SAMPLE}.vcf.gz" | \
        awk -v s="$SAMPLE" '
/^#/{ next }
{
    split($8, info, ";")
    af_str = ""
    for(i in info) {
        if(info[i] ~ /^AF=/) {
            af_str = substr(info[i], 4)
            break
        }
    }
    split($5, alts, ",")
    split(af_str, afs, ",")
    for(i = 1; i <= length(alts); i++) {
        af = afs[i]
        if(af == "") af = "NA"
        print s "\t" $1 "\t" $2 "\t" $4 "\t" alts[i] "\t" af
    }
}' >> "$COLLAPSED"
    done
fi