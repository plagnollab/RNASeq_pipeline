# Script to process quantseq data 
PIPEBASE=/SAN/vyplab/HuRNASeq/RNASeq_pipeline
#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

inputDir=/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/threeprimeseq/ko/
support=/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/threeprimeseq/ko/ko_support.tab 

randomPrimingScript=${PIPEBASE}/polyA_random_priming.py
thinningScript=${PIPEBASE}/polyA_remove_weakest.py

outFolder=/SAN/vyplab/IoN_RNAseq/Nicol/ko/
clusterList=${outFolder}/cluster_list.tab

report=${outFolder}/ko_quantseq_report.tab
species=mouse



minReads=5 # minimum number of supporting reads a cluster can have in one sample
mapQual=10 # for  Samtools - what is the probability of the alignment not being the best possible p_wrong= 10^(-mapq/10). Therefore mapQ=10 is 10% chance
minSamples=2 # minimum number of samples a cluster can be found in to be kept
downstreamFlank=1000 # number of nt downstream of each gene region to add for annotation
minProportion=0.05 # minimum contribution a cluster can make to the total

if [ ! -e $outFolder ];then
   mkdir $outFolder
fi

if [ "$species" == "mouse" ];then
   genome=/SAN/vyplab/HuRNASeq/reference_datasets/mm10.fa
fi

# overwrite if exists
# > $clusterList



# echo "QuantSeq pipeline" > $report

# # Convert bam file into bed and filter out multimapped reads 

# i=0
# while read line ; 
# do 
#    set $line 
#    sample=$1 
#    if [ $sample == "sample" ]; then 
#        continue 
#    fi 
   
#    bamfile=${inputDir}/${sample}/${sample}_unique.bam  
   
#    outDir=${outFolder}/${sample}
#    if [ ! -e $outDir ];then
#       mkdir $outDir
#    fi

#    bedfile=${outDir}/${sample}_clusters.bed
#    tmpfile=${outDir}/${sample}_clean.bam

#    echo `date` >> $report
#    echo $sample >> $report
#    echo "========" >> $report
#    echo "Total reads:" `samtools view -c $bamfile` >> $report
   
#    bedfiles[$i]=$bedfile 
#    i=`expr $i + 1`
#    echo $bamfile 
#    echo $bedfile  
#    # q10 is samtools command for a certain quality that includes multimapping - double check this!
#    samtools view -h -q $mapQual $bamfile | 
# 	  awk '$6 !~ /N/ || $1 ~ /^@/ || $1 ~ /^ID/' | 
# 	  sed 's/^ID/@ID/' | 
# 	  samtools view -bh - > $tmpfile


#    echo "Uniquely mapped:" `samtools view -c $tmpfile` >> $report 
   	
# 	# make bed and count occurences - filter out anything with fewer than 5 reads - possible to change
# 	bedtools bamtobed -i $tmpfile | 
#       bedtools merge -s -c 4,3,6 -o count,count_distinct,distinct | 
#       awk -v minReads=$minReads '$4 >= minReads { print $0 }' -  > $bedfile

#    echo $bedfile >> $clusterList

#    echo "Bed clusters:\t" `wc -l $bedfile | awk '{print $1}'` >> $report

#    	#rm $tmpfile
# done < $support

# FEATURES TO ADD
# 	remove random priming artifacts (high A content 10nt downstream but not canonical polyA sequence)
#	remove polyAs that contribute <5% of total counts in a gene

# how many polyAs are within or close to annotated UTRs? check 


# Combine all the individual bed files into one merged bed

# bedlist=$(printf " %s" "${bedfiles[@]}")
# # collapse bash array
# bedlist=${bedlist:1}

# echo $bedlist 
allBed=${outFolder}/all_samples.bed
allBedAnno=${outFolder}/all_samples_annotated.bed
# # check if in at least 2 samples
cat $clusterList | 
   xargs cat |
#cat $bedlist | 
   bedtools sort | 
   bedtools merge -s -c 4,3,6 -o sum,count,distinct | 
   awk -v minSamples=$minSamples '$5 >= minSamples' > ${allBed}

echo "total merged clusters: " `wc -l ${allBed} | awk '{print $1}' ` #>> $report



echo "annotating the bed file" 

# maybe just use the gene coordinates?
#annotation=/SAN/vyplab/HuRNASeq/reference_datasets/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gtf
annotation=/SAN/vyplab/HuRNASeq/reference_datasets/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed_genes_only.gtf

# # TODO: flank every gene coord downstream to catch unannotated polyAs
# bedtools intersect  -wb -a $annotation -b $allBed > ${allBedAnno}
# echo "clusters that intersect a gene body: " `wc -l $allBedAnno | awk '{print $1}' ` >> $report

# bedtools intersect -S -wb -a $annotation -b $allBed > ${allBedAnno}_stranded
# echo "clusters that intersect a gene body of the correct strand: " `wc -l ${allBedAnno}_stranded | awk '{print $1}' ` >> $report


 # instead intersect with a 1000bp window downstream of the gene body
bedtools window -sw -l 0 -r $downstreamFlank -Sm -a $annotation -b $allBed > ${allBedAnno}_stranded_${downstreamFlank}_window 
echo "clusters that intersect a gene body of the correct strand with ${downstreamFlank} bp downstream window: " `wc -l ${allBedAnno}_stranded_${downstreamFlank}_window | awk '{print $1}' ` #>> $report

# find ensembl gene name from gtf intersect
awk -F'"' '{print $2}' $allBedAnno > tmp.tab
# sort and remove duplicate rows 
cut -f10-15 ${allBedAnno}_stranded_${downstreamFlank}_window | 
   paste - tmp.tab | 
   sort -k1,1 -k2,2n -u - |
   awk 'NF == 7'  > tmp2.bed # sanity check 

#bedtools sort -i tmp2.bed  > tmp3.bed 
mv tmp2.bed ${allBedAnno}_stranded_${downstreamFlank}_window.bed
rm tmp.tab 


#exit

echo "removing random priming artifacts" 

 cleanBed=${allBedAnno}_stranded_${downstreamFlank}_window_clean_new_method.bed

python $randomPrimingScript \
   ${allBedAnno}_stranded_${downstreamFlank}_window.bed \
   $genome \
   $cleanBed

echo "After removal of random priming artifacts: " `wc -l $cleanBed | awk '{print $1}'` #>> $report


#exit

echo "thinning out weakest clusters"

python $thinningScript \
   $cleanBed \
   $minProportion \
   ${cleanBed}_thinned

echo "After removal of weakest clusters (min proportion >=" $minProportion "):" `wc -l ${cleanBed}_thinned | awk '{print $1}'`

# Note: the resulting annotated bed files will have more lines than the original clusters bed file 
# Each cluster location can be annotated as more than one gene 

# each sample has a tmp.bam - keep for the next step
exit

# Now use the merge bed file to count 
while read line ; 
do 
   set $line 
   sample=$1 
   if [ $sample == "sample" ]; then 
       continue 
   fi 
   
   bamfile=${inputDir}/${sample}/${sample}_unique.bam  
   #bedfile=${inputDir}/${sample}/${sample}_clusters.bed
   obedfile=${outFolder}/${sample}/${sample}_clusters_counts.bed 
   tmpfile=${outFolder}/${sample}/${sample}_clean.bam 
   echo $obedfile  
   #echo $bedfile 
   echo $allBed
   #samtools view -h -q10 $bamfile | awk '$6 !~ /N/ || $1 ~ /^@/ || $1 ~ /^ID/' | sed 's/^ID/@ID/' | samtools view -bh - > $tmpfile 
   bedtools intersect -c -wa -a $cleanBed -b $tmpfile -c > $obedfile 
   #rm $tmpfile 
done < $support 
