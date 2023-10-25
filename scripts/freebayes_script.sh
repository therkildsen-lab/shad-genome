#!/bin/bash

BAMLIST=$1
REF=$2
THREADS=$3
OUTDIR=$4
POPNAME=$5clear

freebayes-parallel <(fasta_generate_regions.py $REF 100000) $THREADS -f $REF -L $BAMLIST | vcffilter -f "QUAL > 20" | bgzip > $OUTDIR'/'$POPNAME'_freebayes_results.vcf.gz'

#freebayes -f /workdir/azwad/shad_genome/genomes/fAloSap1/assembly/Alosa_sapidissima_fAloSap1.pri_genome.fna fAloSap1_mapped.sorted.bam |
#vcffilter -f "QUAL > 20" > results.vcf
