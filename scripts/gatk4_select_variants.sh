#!/bin/bash

INPUT_VCF=$1
SAMPLE_NAME=$(basename "$INPUT_VCF" .vcf.gz)

singularity exec --bind $PWD --pwd $PWD /programs/gatk-4.5.0.0/gatk.sif \
gatk SelectVariants \
-V $INPUT_VCF \
-select-type SNP \
-O "$SAMPLE_NAME"_snp.vcf.gz

singularity exec --bind $PWD --pwd $PWD /programs/gatk-4.5.0.0/gatk.sif \
gatk SelectVariants \
-V $INPUT_VCF \
-select-type INDEL \
-O "$SAMPLE_NAME"_indel.vcf.gz

singularity exec --bind $PWD --pwd $PWD /programs/gatk-4.5.0.0/gatk.sif \
gatk SelectVariants \
-V $INPUT_VCF \
-select-type NO_VARIATION \
-O "$SAMPLE_NAME"_invariant.vcf.gz
