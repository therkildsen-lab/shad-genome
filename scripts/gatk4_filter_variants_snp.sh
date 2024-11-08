#!/bin/bash

INPUT_VCF=$1
SAMPLE_NAME=$(basename "$INPUT_VCF" .vcf.gz)

singularity exec --bind $PWD --pwd $PWD /programs/gatk-4.5.0.0/gatk.sif \
gatk VariantFiltration \
-V $INPUT_VCF \
-filter "DP > 83" --filter-name "DP83" \
-filter "DP < 13" --filter-name "DP13" \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 60" --filter-name "QUAL60" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O "$SAMPLE_NAME"_filtered.vcf.gz