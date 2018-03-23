#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# extract polyA clusters
echo "writing cleaned BAM files"

#sh $createCleanedBams --support ${support} --outFolder ${outFolder} --code ${code} --clusterList ${clusterList}
until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
   --support)
		shift
		support=$1;;
   --iFolder)
		shift
		iFolder=$1;;
   --outFolder)
     shift
     outFolder=$1;;
   --code)
        shift
        code=$1;;
   --clusterList)
       shift
        clusterList=$1;;
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
   

  outDir=${outFolder}/cleaned_bams/${code}

  if [ ! -e $outDir ];then
      mkdir -p $outDir
  fi

  bamfile=${outFolder}/$sample/${sample}_clean.bam
  cleanedBam=${outDir}/${sample}_${code}.bam

  if [ ! -e $bamfile ];then
      echo $bamfile does not exist
      exit 1
  fi

  # check that cluster file given exists
  if [ ! -e $clusterList ];then
      echo $clusterList does not exist
      exit 1
  fi

  # extract all the reads in a bam file that overlap a given list of clusters  
   #
   echo writing cleaned bam to $cleanedBam
   samtools view -bh $bamfile -L $clusterList > $cleanedBam
   echo indexing
   samtools index $cleanedBam

done < $support