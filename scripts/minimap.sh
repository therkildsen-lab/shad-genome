
#!/bin/bash
# minimap script

#1:sr/map-pb/map-ont/map-hifi
#2:path/to/assembly
#3:output name
#4:input fastq files

export PATH=/programs/minimap2-2.26:$PATH

minimap2 -ax $1 -t 16 -R @RG\\tID:1\\tPL:PACBIO\tLB:LB1\\tSM:fAloSap1 $2 ${@:4} | samtools sort -@16 -O BAM -o $3.bam 
