#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# count clusters

until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
   --masterList)
        shift
        masterList=$1;;
   --support)
       shift
        support=$1;;
   --outFolder)
       shift
        outFolder=$1;;
   --iFolder)
		shift
		iFolder=$1;;
   -* )
       stop "Unrecognized option: $1"
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
    echo $1
done

# masterList : contains all the clusters we want to count
# support : contains list of sample IDs
# outFolder : the same outFolder that has been used for all the previous steps

# check that all the cleaned cluster bams still exist
clusterBams=`awk -v outFolder=${outFolder} 'NR>1{print outFolder"/"$1"/"$1"_clean.bam"}' $support ` 

fail=0
for bam in $clusterBams; do
	if [ ! -e $bam ];then
		echo $bam cannot be found
		fail=1
	else
		echo $bam found
	fi
done

if [ fail == 1 ]; then
	exit 1
fi

# 
echo "counting clusters against master list in" $masterList

# Now use the merge bed file to count 
while read line ; 
do 
   set $line 
   sample=$1 
   if [ $sample == "sample" ]; then 
       continue 
   fi 
   
   bamfile=${iFolder}/${sample}/${sample}_unique.bam  
   #bedfile=${iFolder}/${sample}/${sample}_clusters.bed
   obedfile=${outFolder}/${sample}/${sample}_clusters_counts.bed 
   tmpfile=${outFolder}/${sample}/${sample}_clean.bam 
   #echo $obedfile  
   #echo $bedfile 

   #samtools view -h -q10 $bamfile | awk '$6 !~ /N/ || $1 ~ /^@/ || $1 ~ /^ID/' | sed 's/^ID/@ID/' | samtools view -bh - > $tmpfile 
   bedtools intersect -c -wa -a $masterList -b $tmpfile -c > $obedfile
   echo "writing " $obedfile 
   #rm $tmpfile 
done < $support 