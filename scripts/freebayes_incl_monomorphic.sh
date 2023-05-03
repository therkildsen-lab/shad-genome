#!/bin/bash

INPUT_BAM=$1
THREADS=$2

freebayes-parallel <(fasta_generate_regions.py ../assembly/Alosa_sapidissima_fAloSap1.pri_genome.fna.fai 100000) $THREADS -f ../assembly/Alosa_sapidissima_fAloSap1.pri_genome.fna --report-monomorphic $INPUT_BAM | vcffilter -f "QUAL > 20" > freebayes_results_incl_monomorphic.vcf

gzip freebayes_results_incl_monomorphic.vcf
