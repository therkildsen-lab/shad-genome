#!/bin/bash

INPUT_BAM=$1
THREADS=$2

freebayes-parallel <(fasta_generate_regions.py ../assembly/Alosa_sapidissima_fAloSap1.pri_genome.fna.fai 100000) $THREADS -f ../assembly/Alosa_sapidissima_fAloSap1.pri_genome.fna $INPUT_BAM | vcffilter -f "QUAL > 20" > freebayes_results.vcf

gzip freebayes_results.vcf

#freebayes -f /workdir/azwad/shad_genome/genomes/fAloSap1/assembly/Alosa_sapidissima_fAloSap1.pri_genome.fna fAloSap1_mapped.sorted.bam | 
#vcffilter -f "QUAL > 20" > results.vcf
