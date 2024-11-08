#!/bin/bash

REF=$1
INPUT_VCF=$2
SAMPLE_NAME=$(basename "$INPUT_VCF" .vcf.gz)

singularity exec --bind $PWD --pwd $PWD /programs/gatk-4.5.0.0/gatk.sif \
gatk GenotypeGVCFs \
-R $REF \
-all-sites \
-V $INPUT_VCF \
-O "$SAMPLE_NAME"_genotypegvcfs_unfiltered.vcf.gz
