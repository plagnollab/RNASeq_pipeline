#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# sh $mergeClusters --outFolder ${outFolder} \
#                   --code ${code} \
#                   --clusterList ${clusterList} \
#                   --minSamples ${minSamples} \
#                   --annotation ${annotation} \
#                   --downstreamFlank ${downstreamFlank}

# merge clusters together
echo "merging clusters"


until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
   --minSamples)
		shift
		minSamples=$1;;
   --annotation)
		shift
		annotation=$1;;
   --downstreamFlank)
		shift
		downstreamFlank=$1;;
   --outFolder)
        shift
        outFolder=$1;;
   --code)
        shift
        code=$1;;
   --clusterList)
       shift
        clusterList=$1;;
   --mergedOut)
		shift
		mergedOut=$1;;
   -* )
       stop "Unrecognized option: $1"
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
    echo $1
done



# echo $bedlist 
#allBed=${outFolder}/all_samples.bed
#allBedAnno=${outFolder}/all_samples_annotated.bed
# # check if in at least 2 samples
cat $clusterList | 
	   xargs cat |
   bedtools sort | 
   bedtools merge -s -c 4,3,6 -o sum,count,distinct | 
   awk -v minSamples=$minSamples '$5 >= minSamples' > ${mergedOut}_raw

echo "total merged clusters:	"`wc -l ${mergedOut}_raw | awk '{print $1}' ` #>> $report


echo "annotating the bed file" 

# maybe just use the gene coordinates?
#annotation=/SAN/vyplab/HuRNASeq/reference_datasets/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gtf

# # TODO: flank every gene coord downstream to catch unannotated polyAs
# bedtools intersect  -wb -a $annotation -b $allBed > ${allBedAnno}
# echo "clusters that intersect a gene body: " `wc -l $allBedAnno | awk '{print $1}' ` >> $report

# bedtools intersect -S -wb -a $annotation -b $allBed > ${allBedAnno}_stranded
# echo "clusters that intersect a gene body of the correct strand: " `wc -l ${allBedAnno}_stranded | awk '{print $1}' ` >> $report


 # instead intersect with a 1000bp window downstream of the gene body
bedtools window -sw -l 0 -r ${downstreamFlank} -Sm -a ${annotation} -b ${mergedOut}_raw > ${mergedOut}_annotated
# find ensembl gene name from gtf intersect
awk -F'"' '{print $2}' ${mergedOut}_annotated > ${mergedOut}_geneIDs

# sort and remove duplicate rows 
cut -f10-15 ${mergedOut}_annotated | 
   	  paste - ${mergedOut}_geneIDs | 
      		sort -k1,1 -k2,2n -u - |
   			awk 'NF == 7'  > ${mergedOut}.bed # sanity check 


echo "clusters that intersect a gene body of the correct strand with ${downstreamFlank} bp downstream window: " `wc -l ${mergedOut}.bed | awk '{print $1}' ` #>> $report

rm ${mergedOut}_raw ${mergedOut}_annotated ${mergedOut}_geneIDs 



