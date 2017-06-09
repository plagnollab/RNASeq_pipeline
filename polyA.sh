# Script to process quantseq data 

# #!/bin/bash

inputDir=/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/threeprimeseq/ko/
support=/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/threeprimeseq/ko/ko_support.tab 

# Convert bam file into bam and filter out multimapped reads 

i=0
while read line ; 
do 
   set $line 
   sample=$1 
   if [ $sample == "sample" ]; then 
       continue 
   fi 
   
   bamfile=${inputDir}/${sample}/${sample}_unique.bam  
   bedfile=${inputDir}/${sample}/${sample}_clusters.bed
   tmpfile=${inputDir}/${sample}/tmp.bam 

   bedfiles[$i]=$bedfile 
   i=`expr $i + 1`
   echo $bamfile 
   echo $bedfile  
   #samtools view -h -q10 $bamfile | awk '$6 !~ /N/ || $1 ~ /^@/ || $1 ~ /^ID/' | sed 's/^ID/@ID/' | samtools view -bh - > $tmpfile 
   #bedtools bamtobed -i $tmpfile | bedtools merge -s -c 4,3,6 -o count,count_distinct,distinct | awk '$4 >= 5 { print $0 }' -  > $bedfile
   #rm $tmpfile
done < $support

# Now combine all the individual bed files into one 
bedlist=$(printf " %s" "${bedfiles[@]}")
bedlist=${bedlist:1}

echo $bedlist 
allBed=${inputDir}/all_samples.bed
allBedAnno=${inputDir}/all_samples_annotated.bed
cat $bedlist | bedtools sort | bedtools merge -s -c 4,3,6 -o sum,count,distinct | awk '$5 > 2' > ${inputDir}/all_samples.bed

echo "annotating the bed file" 
annotation=/SAN/vyplab/HuRNASeq/reference_datasets/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gtf
bedtools intersect  -wb -a $annotation -b $allBed > $allBedAnno
awk -F'"' '{print $2}' $allBedAnno > tmp.tab 
cut -f10-15 $allBedAnno | paste - tmp.tab | sort -u - | bedtools sort  > tmp2.bed 
mv tmp2.bed $allBedAnno
rm tmp.tab 

# Note: the resulting annotated bed files will have more lines than the original clusters bed file 
# Each cluster location can be annotated as more than one gene 

# Now use the merge bed file to count 
while read line ; 
do 
   set $line 
   sample=$1 
   if [ $sample == "sample" ]; then 
       continue 
   fi 
   
   bamfile=${inputDir}/${sample}/${sample}_unique.bam  
   bedfile=${inputDir}/${sample}/${sample}_clusters.bed
   obedfile=${inputDir}/${sample}/${sample}_clusters_counts.bed 
   tmpfile=${inputDir}/${sample}/tmp.bam 
   echo $obedfile  
   echo $bedfile 
   echo $allBed
   #samtools view -h -q10 $bamfile | awk '$6 !~ /N/ || $1 ~ /^@/ || $1 ~ /^ID/' | sed 's/^ID/@ID/' | samtools view -bh - > $tmpfile 
   bedtools intersect -c -wa -a $allBedAnno -b $tmpfile -c > $obedfile 
   rm $tmpfile 
done < $support 
