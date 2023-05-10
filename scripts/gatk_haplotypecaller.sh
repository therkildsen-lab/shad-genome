#!/bin/bash

REF=$1
INPUT_BAM=$2

java -jar /programs/bin/GATK/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -nct 16 \
    -R $REF \
    -I $INPUT_BAM \
    -ERC BP_RESOLUTION \
    -out_mode EMIT_ALL_SITES \
    --genotyping_mode DISCOVERY \
    --allow_potentially_misencoded_quality_scores \
    -o reads/gatk_raw_variants.vcf
