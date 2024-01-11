#!/bin/bash

REF=$1
BAM=$2
MINDEPTH=$3
MAXDEPTH=$4

filename=$(basename -- "$BAM")
filename="${filename%.*}"

bcftools mpileup -Ou -X ont --threads 16 -f $REF $BAM | bcftools call -P 0.01 -c - | vcfutils.pl vcf2fq -d $MINDEPTH -D $MAXDEPTH | gzip > "$filename"_diploid.fq.gz
