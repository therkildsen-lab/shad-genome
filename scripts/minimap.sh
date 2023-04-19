
#!/bin/bash
# minimap script

#1:sr/map-pb/map-ont
#2:path/to/assembly
#3:output name
#4:input fastq files

minimap2 -ax $1 -t 8 $2 ${@:4} > $3.sam

samtools sort -@8 -O BAM -o $3.bam  $3.sam
