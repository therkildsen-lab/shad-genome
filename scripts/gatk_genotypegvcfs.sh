#!/bin/bash

REF=$1

export JAVA_HOME=/usr/local/jdk1.8.0_121
export PATH=$JAVA_HOME/bin:$PATH
module load java/1.8.0

java -Xmx100G -jar /programs/bin/GATK/GenomeAnalysisTK.jar \
    -T GenotypeGVCFs \
    -R $REF \
    -allSites \
    -stand_call_conf 0 \
    --variant reads/gatk_raw_variants.vcf \
    -o reads/gatk_processed_variants.vcf

