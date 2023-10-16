#!/bin/bash

VCF=$1
MINDEPTH=$2
MAXDEPTH=$3

filename=$(basename -- "$VCF")
filename="${filename%.*}"

vcfutils.pl vcf2fq -d $MINDEPTH -D $MAXDEPTH $VCF | gzip > "$filename"_diploid.fq.gz