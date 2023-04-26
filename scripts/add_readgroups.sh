java -jar /programs/picard-tools-2.26.1/picard.jar AddOrReplaceReadGroups \
       I=$1 \
       O=$2 \
       RGID=1 \
       RGLB=lib1 \
       RGPL=PACBIO \
       RGPU=unit1 \
       RGSM=fAloSap1