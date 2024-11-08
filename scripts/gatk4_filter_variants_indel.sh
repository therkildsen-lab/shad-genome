#!/bin/bash

INPUT_VCF=$1
SAMPLE_NAME=$(basename "$INPUT_VCF" .vcf.gz)

singularity exec --bind $PWD --pwd $PWD /programs/gatk-4.5.0.0/gatk.sif \
gatk VariantFiltration \
-V $INPUT_VCF \
-filter "DP > 83" --filter-name "DP83" \
-filter "DP < 13" --filter-name "DP13" \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 60.0" --filter-name "QUAL60" \
-filter "FS > 200.0" --filter-name "FS200" \
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
-O "$SAMPLE_NAME"_filtered.vcf.gz