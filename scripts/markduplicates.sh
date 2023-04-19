BAM=$1
JAVA=${3:-java} # Path to java
PICARD=${4:-/programs/picard-tools-2.9.0/picard.jar} # Path to picard tools
BAMUTIL=${5:-/programs/bamUtil/bam} # Path to bamUtil

SAMPLEPREFIX=`echo ${BAM%.bam}`

## Remove duplicates and print dupstat file
	# We used to be able to just specify picard.jar on the CBSU server, but now we need to specify the path and version
	$JAVA -Xmx60g -jar $PICARD MarkDuplicates I=$BAM O=$SAMPLEPREFIX'_dedup.bam' M=$SAMPLEPREFIX'_dupstat.txt' VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
	
