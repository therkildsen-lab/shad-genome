#!/bin/bash

export PATH=/programs/psmc-0.6.5:$PATH

DIPLOID_FASTQ=$1
filename=$(basename -- "$DIPLOID_FASTQ")
filename="${filename%.*.*}"

/workdir/azwad/shad-genome/psmc/psmc/utils/fq2psmcfa -q20 $DIPLOID_FASTQ > $filename.psmcfa

/workdir/azwad/shad-genome/psmc/psmc/utils/splitfa $filename.psmcfa > $filename'_split.psmcfa'

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o $filename.psmc $filename.psmcfa

seq 100 | xargs -i echo psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" \
-o round-{}.psmc $filename'_split.psmcfa' | sh

cat $filename.psmc round-*.psmc > $filename'_combined.psmc'

rm round-*.psmc