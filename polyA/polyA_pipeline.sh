#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Pipeline to process quantseq data 
# PIPEBASE=/SAN/vyplab/HuRNASeq/RNASeq_pipeline/polyA

# # for testing
# iFolder=/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/threeprimeseq/ko/
# #support=/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/threeprimeseq/ko/ko_support.tab 
# support=/SAN/vyplab/IoN_RNAseq/Nicol/ko/test_support.tab
# outFolder=/SAN/vyplab/IoN_RNAseq/Nicol/ko/
# code="FUS_ko"
# species=mouse
# submit=no

# case statement to get variables

until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
   --species)
       shift
       species=$1;;
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
   --submit)
       shift
       submit=$1;;
   --PIPEBASE)
       shift
       PIPEBASE=$1;;
   -* )
       stop "Unrecognized option: $1"
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
    echo $1
done


#languages
R=/share/apps/R-3.3.2/bin/R
python=/share/apps/python-2.7.8/bin/python
export PYTHONPATH=/SAN/vyplab/HuRNASeq:/SAN/vyplab/HuRNASeq/python-2.7.8/lib/python2.7/site-packages

# scripts
extractPolyAClusters=${PIPEBASE}/extractPolyAClusters.sh
mergeClusters=${PIPEBASE}/mergeClusters.sh
randomPrimingScript=${PIPEBASE}/polyA_random_priming.py
thinningScript=${PIPEBASE}/polyA_remove_weakest.py
countClusters=${PIPEBASE}/countClusters.sh
createCleanedBams=${PIPEBASE}/createCleanedBams.sh
polyA_dexseq=${PIPEBASE}/polyA_dexseq.R
analyseClusters=${PIPEBASE}/analyseClusters.R

# variables
minReads=5 # minimum number of supporting reads a cluster can have in one sample
mapQual=10 # for  Samtools - what is the probability of the alignment not being the best possible p_wrong= 10^(-mapq/10). Therefore mapQ=10 is 10% chance
minSamples=2 # minimum number of samples a cluster can be found in to be kept
downstreamFlank=1000 # number of nt downstream of each gene region to add for annotation
minProportion=0.05 # minimum contribution a cluster can make to the total
nCores=8 # number of cores for DEXSeq
fail=0

if [ "$species" == "mouse" ];then
   genome=/SAN/vyplab/HuRNASeq/reference_datasets/mm10.fa
   annotation=/SAN/vyplab/HuRNASeq/reference_datasets/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed_genes_only.gtf
   biomartAnnotation=/SAN/vyplab/HuRNASeq/reference_datasets/RNASeq/Mouse/biomart_annotations_mouse.tab
elif [ "$species" == "human" ];then
   genome=/SAN/vyplab/HuRNASeq/reference_datasets/hg38.fa
   annotation=/SAN/vyplab/HuRNASeq/reference_datasets/RNASeq/Human_hg38/Homo_sapiens.GRCh38.82_fixed.gtf
   biomartAnnotation=/SAN/vyplab/HuRNASeq/reference_datasets/RNASeq/Human_hg38/biomart_annotations_human.tab
else
   echo $species is not a valid species
   fail=1
fi

# make folders
for folder in ${outFolder}/cluster/submission ${outFolder}/cluster/out ${outFolder}/cluster/error;do
   if [ ! -e $folder ];then
      mkdir -p $folder
   fi
done

# check everything exists

for file in $extractPolyAClusters $mergeClusters $randomPrimingScript $thinningScript $countClusters $outFolder $genome $annotation;do
   if [ ! -e ${file} ]; then
      echo $file not found!
      fail=1
   fi
done

# check bam files exist
bamFiles=`awk -v iFolder=${iFolder} 'NR>1{print iFolder"/"$1"/"$1"_unique.bam"}' $support` 
for bam in $bamFiles; do
   if [ ! -e $bam ]; then
      echo ${bam} cannot be found
      fail=1
   else
      echo ${bam} found
   fi
done

if [ $fail == 1 ];then
   exit 0
fi




jobScript=${outFolder}/cluster/submission/polyA_submission.sh
dexseqScript=${outFolder}/cluster/submission/polyA_dexseq.sh

clusterList=${outFolder}/${code}_cluster_list.tab
mergedOut=${outFolder}/${code}_all_samples

echo writing pipeline script to $jobScript

echo "
#$ -S /bin/bash
#$ -l h_vmem=8G,tmem=8G
#$ -l h_rt=72:00:00
#$ -pe smp 1
#$ -R y
#$ -o ${outFolder}/cluster/out
#$ -e ${outFolder}/cluster/error
#$ -N polyA_pipeline_${code}
#$ -wd ${outFolder}
echo \$HOSTNAME >&2

# overwrite if exists; if not create empty file
 > $clusterList

sh $extractPolyAClusters --support ${support} \\
                         --outFolder ${outFolder} \\
                         --iFolder ${iFolder} \\
                         --code ${code} \\
                         --clusterList ${clusterList} \\
                         --minReads ${minReads} \\
                         --mapQual ${mapQual} 

# # Convert bam file into bed and filter out multimapped reads 


sh $mergeClusters --outFolder ${outFolder} \\
                  --code ${code} \\
                  --clusterList ${clusterList} \\
                  --minSamples ${minSamples} \\
                  --annotation ${annotation} \\
                  --downstreamFlank ${downstreamFlank} \\
                  --mergedOut ${mergedOut}

# remove clusters that probably result from random priming

$python $randomPrimingScript ${mergedOut}.bed \\
                             $genome \\
                             ${mergedOut}_artifacts_removed.bed

echo \"After removal of random priming artifacts: \" \`wc -l ${mergedOut}_artifacts_removed.bed | awk '{print \$1}'\` 

# create bam files with random priming removed
sh $createCleanedBams --support ${support} \\
                      --outFolder ${outFolder} \\
                      --iFolder ${iFolder} \\
                      --code random_priming_only \\
                      --clusterList ${mergedOut}_artifacts_removed.bed \\


$python $thinningScript \\
   ${mergedOut}_artifacts_removed.bed \\
   $minProportion \\
   ${mergedOut}_artifacts_removed_thinned.bed

echo \"After removal of weakest clusters (min proportion >=\" $minProportion \"):\" \`wc -l ${mergedOut}_artifacts_removed_thinned.bed | awk '{print \$1}'\`

# create bam files with random priming removed AND only the most contributing clusters
sh $createCleanedBams --support ${support} \\
                      --outFolder ${outFolder} \\
                      --iFolder ${iFolder} \\
                      --code most_important_clusters_${minProportion} \\
                      --clusterList ${mergedOut}_artifacts_removed_thinned.bed \\


# Note: the resulting annotated bed files will have more lines than the original clusters bed file 
# Each cluster location can be annotated as more than one gene 

# each sample has a tmp.bam - keep for the next step

# count with the unthinned and the thinned annotation

sh $countClusters --support ${support} \\
                  --outFolder ${outFolder} \\
                  --iFolder ${iFolder} \\
                  --masterList ${mergedOut}_artifacts_removed.bed \\
                  --mode no_artifacts

sh $countClusters --support ${support} \\
                  --outFolder ${outFolder} \\
                  --iFolder ${iFolder} \\
                  --masterList ${mergedOut}_artifacts_removed_thinned.bed \\
                  --mode thinned
" > $jobScript

echo writing DEXSeq code to $dexseqScript

# run dexseq on thinned cluster list
echo "
#$ -S /bin/bash
#$ -l h_vmem=5G,tmem=5G
#$ -l h_rt=72:00:00
#$ -pe smp ${nCores}
#$ -R y
#$ -o ${outFolder}/cluster/out
#$ -e ${outFolder}/cluster/error
#$ -N polyA_dexseq_${code}
#$ -wd ${outFolder}
echo \$HOSTNAME >&2


${R}script ${polyA_dexseq} --support.tab $support \\
                           --code $code \\
                           --output.dir ${outFolder} \\
                           --input.dir $iFolder \\
                           --biomartAnnotation $biomartAnnotation \\
                           --nCores $nCores \\
                           --mode no_artifacts


${R}script ${polyA_dexseq} --support.tab $support \\
                           --code $code \\
                           --output.dir ${outFolder} \\
                           --input.dir $iFolder \\
                           --biomartAnnotation $biomartAnnotation \\
                           --nCores $nCores \\
                           --mode thinned


# do cluster analysis

${R}script ${analyseClusters} --code $code \\
                              --output.dir ${outFolder}/results \\
                              --biomartAnnotation $biomartAnnotation \\
                              --mode no_artifacts \\
                              --dexseqResults ${outFolder}/results/dexseq_results_${code}_no_artifacts.tab

${R}script ${analyseClusters} --code $code \\
                              --output.dir ${outFolder}/results \\
                              --biomartAnnotation $biomartAnnotation \\
                              --mode thinned \\
                              --dexseqResults ${outFolder}/results/dexseq_results_${code}_thinned.tab

" > $dexseqScript

if [ "$submit" == "yes" ];then
   qsub $jobScript
   qsub -hold_jid polyA_pipeline_${code} $dexseqScript
fi

