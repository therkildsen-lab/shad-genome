#!/bin/bash

REF=$1
BAM=$2
filename=$(basename -- "$BAM")
filename="${filename%.*}"

bcftools mpileup -Ou -X pacbio-ccs --threads 16 -f $REF $BAM | bcftools call --threads 16 -c -Ov -o $filename"_for_diploid.vcf"
