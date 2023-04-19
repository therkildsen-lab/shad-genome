#!/bin/bash

REF=$1
BAM=$2
filename=$(basename -- "$BAM")
filename="${filename%.*}"

bcftools mpileup -C50 -Ou -f $REF $BAM | bcftools call -c - | vcfutils.pl vcf2fq -d 13 -D 82 | gzip > "$filename"_diploid.fq.gz
