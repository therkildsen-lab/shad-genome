#!/bin/bash

REF=$1
INPUT_BAM=$2
SAMPLE_NAME=$(basename "$INPUT_BAM" .bam)

singularity exec --bind $PWD --pwd $PWD /programs/gatk-4.5.0.0/gatk.sif gatk HaplotypeCaller \
-R $REF \
-I $INPUT_BAM \
-O "$SAMPLE_NAME"_haplotypecaller.vcf.gz \
-ERC BP_RESOLUTION \
-G StandardAnnotation \
-G AS_StandardAnnotation