#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# extract polyA clusters
echo "extracting polyA clusters"

#sh $extractPolyAClusters --support ${support} --outFolder ${outFolder} --code ${code} --clusterList ${clusterList}
until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
   --support)
		shift
		support=$1;;
   --outFolder)
        shift
        outFolder=$1;;
   --iFolder)
		shift
		iFolder=$1;;
   --code)
        shift
        code=$1;;
   --clusterList)
       shift
        clusterList=$1;;
   --minReads)
		shift
		minReads=$1;;
	--mapQual)
		shift
		mapQual=$1;;
   -* )
       stop "Unrecognized option: $1"
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
    echo $1
done


i=0
while read line ; 
do 
   set $line 
   sample=$1 
   if [ $sample == "sample" ]; then 
       continue 
   fi 

   echo $sample
   
   bamfile=${iFolder}/${sample}/${sample}_unique.bam 

   if [ ! -e $bamfile ];then
   		echo $bamfile does not exist
   		exit 1
   fi

   outDir=${outFolder}/${sample}
   
   if [ ! -e $outDir ];then
      mkdir $outDir
   fi

   bedfile=${outDir}/${sample}_clusters.bed
   tmpfile=${outDir}/${sample}_clean.bam

   echo "	input bam:	"$bamfile 
   echo "	total reads:	"`samtools view -c $bamfile` 
   
   #bedfiles[$i]=$bedfile 
   #i=`expr $i + 1`


   
   # q10 is samtools command for a certain quality that includes multimapping - double check this!
   samtools view -h -q $mapQual $bamfile | 
	  awk '$6 !~ /N/ || $1 ~ /^@/ || $1 ~ /^ID/' | 
	  sed 's/^ID/@ID/' | 
	  samtools view -bh - > $tmpfile

   echo "	Uniquely mapped:	"`samtools view -c $tmpfile`
   echo "	output bed:	"$bedfile  
   	
	# make bed and count occurences - filter out anything with fewer than 5 reads - possible to change
	bedtools bamtobed -i $tmpfile | 
      bedtools merge -s -c 4,3,6 -o count,count_distinct,distinct | 
      awk -v minReads=$minReads '$4 >= minReads { print $0 }' -  > $bedfile

   echo $bedfile >> $clusterList

   echo "	Raw clusters:	"`wc -l $bedfile | awk '{print $1}'`

   	#rm $tmpfile
done < $support
