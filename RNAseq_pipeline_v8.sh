#! /bin/bash
# this script has only been tested with bash and may not work with other shells

# script exits if return value of a command is not zero
set -e
# this forces all variables to be defined
set -u
# for debugging prints out every line before executing it
set -x

# prints to stderr in red
function error() { >&2 echo -e "\033[31m$*\033[0m"; }
function stop() { error "$*"; exit 1; }
try() { "$@" || stop "cannot $*"; }
function files_exist() { for file in $*; do if [ ! -e $file ]; then stop "File $file does not exist" ; fi; done }

#computer=vanHeel
computer=
superLong=

if [[ "$computer" == "vanHeel" ]]
then
    software=/data_n2/vplagnol/Software
    pythonbin=/software/additional/epd-7.3.1/bin/python
    Rbin=/data_n2/vplagnol/Software/R-3.0.2/bin/R
    Rscript=/data_n2/vplagnol/Software/R-3.0.2/bin/Rscript
    dexseqCount=/data_n2/vplagnol/Rlibs/installed/DEXSeq/python_scripts/dexseq_count.py
    javaTemp2="/data_n1/vanheel_singlecellgenomics/tmp"
    javaTemp="TMP_DIR=${javaTemp2}"
else
    software=/cluster/project8/vyp/vincent/Software
    pythonbin=/share/apps/python-2.7.1/bin/python2.7
    if [ ! -e $pythonbin ]; then pythonbin=/share/apps/python-2.7.8/bin/python2.7; fi
    ##Rbin=/cluster/project8/vyp/vincent/Software/R-3.1.2/bin/R
    misoRunEvents=/cluster/project8/vyp/vincent/Software/misopy-0.4.9/misopy/run_events_analysis.py
    runMiso=/cluster/project8/vyp/vincent/Software/misopy-0.4.9/misopy/run_miso.py
    javaTemp2="/scratch2/vyp-scratch2/vincent/java_temp"
    if [ ! -e $javaTemp2 ]; then javaTemp2="/cluster/scratch3/vyp-scratch2/vincent/java_temp/"; fi
    javaTemp="TMP_DIR=${javaTemp2}"
    java=/share/apps/jdk1.7.0_45/bin/java
    if [ ! -e $java ]; then java=/share/apps/jdk1.8.0_25/bin/java; fi
    bigFilesBundleFolder=/cluster/scratch3/vyp-scratch2/reference_datasets

    if [ ! -e $bigFilesBundleFolder ];then
        bigFilesBundleFolder=/SAN/vyplab/HuRNASeq/reference_datasets
    fi

    #which R should we be using?
    #optparse available but just run this to install:
    # > install.packages('getopt')
    # wget --no-check-certificate https://cran.r-project.org/src/contrib/optparse_1.3.2.tar.gz
    # $Rbin CMD install optparse_1.3.2.tar.gz
    Rbin=/cluster/project8/vyp/vincent/Software/R-3.2.2/bin/R
    Rscript=/cluster/project8/vyp/vincent/Software/R-3.2.2/bin/Rscript
    # data.table cannot be installed with this version of R  # Jack: this one works for me! The above does not.
    Rbin=/share/apps/R/bin/R
    Rscript=/share/apps/R/bin/Rscript
    # should fix Rsamtools problem for using DEXSeq
    export LD_LIBRARY_PATH=/share/apps/zlib-1.2.8/lib:$LD_LIBRARY_PATH 
fi


#RNASEQPIPBASE=/Users/pontikos/bin/RNASeq_pipeline/
#RNASEQPIPBASE=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline

echo "Base of RNA-Seq pipeline is located here: $RNASEQPIPBASE"

dexseqCount=${RNASEQPIPBASE}/dexseq_count.py
countPrepareR=${RNASEQPIPBASE}/counts_prepare_pipeline.R
files_exist $countPrepareR
dexseqFinalProcessR=${RNASEQPIPBASE}/dexseq_pipeline_v2.R
deseqFinalProcessR=${RNASEQPIPBASE}/deseq2_pipeline.R
pathwayGOAnalysisR=${RNASEQPIPBASE}/pathwayGO_pipeline.R
topGOAnalysisR=${RNASEQPIPBASE}/topGO_pipeline.R
goseqR=${RNASEQPIPBASE}/goseq_script.R
novosort=${software}/novocraft3/novosort
trimgalore=${software}/trim_galore/trim_galore
cutadapt=/share/apps/python-2.7.8/bin/cutadapt
fastqc=/share/apps/genomics/fastqc-0.11.2/fastqc
#for the old cluster
if [ ! -e $cutadapt ];then cutadapt=/share/apps/python-2.7.6/bin/cutadapt;fi

starexec=/cluster/project8/vyp/vincent/Software/STAR-STAR_2.4.2a/bin/Linux_x86_64_static/STAR
samtools=${software}/samtools-1.2/samtools
rseqQCscripts=${software}/RSeQC-2.3.3/scripts

picardDup=${software}/picard-tools-1.100/MarkDuplicates.jar
picardStats=${software}/picard-tools-1.100/BamIndexStats.jar
picardMetrics=${software}/picard-tools-1.100/CalculateHsMetrics.jar
picardReorder=${software}/picard-tools-1.100/ReorderSam.jar

## should sex chromosomes be kept in the differential expression analysis?
## Niko: if sex matches then yes?
keepSex=FALSE  
force=yes
species=mouse
segmentLength=25

mart=ensembl
db=mmusculus_gene_ensembl
summary=no
prepareCounts=no
Rdexseq=no
Rdeseq=no
RpathwayGO=no
RtopGO=no
oFolder=temp ##default output

#mkdir -p creates folders only if they do not already exist
mkdir -p $oFolder

submit=no
# If QC is wanted, the pipeline will send each fastq or pair of fastqs through FastQC. If adapters are present then the offending fastq files will be trimmed with Trim Galore! and these trimmmed fastqs will be the ones aligned with STAR. 
QC=no
starStep1a=no
starStep1b=no
starStep2=no

stranded=no
libstrand=fr-unstranded
keepDups=FALSE

code=""
dataframe=none
iFolder=""
misoindex=NA
oFolder=RNAseq_processed
stem=""

until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
	--fastqFiles)
	    shift
	    i=0
	    for fileloc in $@; do 
		fastqFiles[ $i ]=$fileloc
		((i=i+1))
	    done;;
	--stranded)
	    shift
	    stranded=yes
	    libstrand=$1;;
	--species)
	    shift
	    species=$1;;
	--force)
	    shift
	    force=$1;;
	--starStep1a)
	    shift
	    starStep1a=$1;;
	--starStep1b)
	    shift
	    starStep1b=$1;;
	--starStep2)
	    shift
	    starStep2=$1;;
	--keep.sex)
	    shift
	    keepSex=$1;;
	--misoindex)
	    shift
	    misoindex=$1;;
	--summary)
	    shift
	    summary=$1;;
	--miso)
	    shift
	    miso=$1;;
	--prepareCounts)
	    shift
	    prepareCounts=$1;;
	--Rdexseq)
	    shift
	    Rdexseq=$1;;
	--Rdeseq)
	    shift
	    Rdeseq=$1;;
	--RpathwayGO)
	    shift
	    RpathwayGO=$1;;
	--RtopGO)
	    shift
	    RtopGO=$1;;
	--iFolder)
	    shift
	    iFolder=$1;;
	--oFolder)
        shift
        oFolder=$1;;
	--dataframe)
        shift
        dataframe=$1;;
	--code)
        shift
        code=$1;;
	--submit)
	    shift
        submit=$1;;
    --stem)
        shift
        stem=$1;;
	--keepDups)
	    shift
	    keepDups=TRUE;;
	--step0_QC)
	    shift
	    step0_QC=$1;;
	--trim_galore)
	    shift
	    trim_galore=$1;;
	--goseq)
		shift
		goseq=$1;;
	-* )
	    stop "Unrecognized option: $1"
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
    echo $1
done


echo "Strand information $stranded $libstrand"

########################## estimate the duration of the jobs

countStrand=no
if [[ "$libstrand" == "fr-firststrand" ]]
then
  countStrand=yes
  countStrandReverse=reverse
elif [[ "$libstrand" == "fr-secondstrand" ]]
then
  countStrand=reverse
  countStrandReverse=yes
else
    echo unknown libstrand $libstrand
fi

if [[ "$stem" == "" ]]
then
    stem=$code
fi

if [[ "$superLong" == "yes" ]]
then
    ((nhours=nhours+nhours))
fi

## create the output folders
clusterFolder=${oFolder}/cluster

mkdir -p ${oFolder} ${clusterFolder} ${clusterFolder}/out ${clusterFolder}/error ${clusterFolder}/R  ${clusterFolder}/submission

files_exist $dataframe 

cp "$dataframe" "${oFolder}/"


case "$species" in
    zebrafish)
        refFolder=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Danio_rerio/NCBI/Zv9
        gtfFile=${refFolder}/Annotation/Genes/genes.gtf    
        fasta=${refFolder}/Sequence/WholeGenomeFasta/genome.fa
        IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome
        gffFile=${RNASEQBUNDLE}/zebrafish/GTF/zebrafish_iGenomes_Zv9_with_ensembl.gff
        annotationFile=${RNASEQBUNDLE}/zebrafish/biomart/biomart_annotations_zebrafish.tab
        #geneModel=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/zebraFish_refSeqTable_zv9_nochr.bed
        #geneModelSummaryStats=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/zebraFish_refSeqTable_zv9_chr1.bed
        ;;
    DvH_sc_human)
        refFolder=/data_n1/vanheel_singlecellgenomics/support/Homo_sapiens/NCBI/build37.2
        gtfFile=${refFolder}/Annotation/Genes/genes.gtf    
        fasta=${refFolder}/Sequence/WholeGenomeFasta/genome.fa
        IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/human_b37_with_spikes
        ;;
    human_hg19_UCSC)
	fasta=/SAN/biomed/biomed14/vyp-scratch/temp_reference/hg19_UCSC/hg19_UCSC.fa
	gtfFile=/SAN/biomed/biomed14/vyp-scratch/temp_reference/hg19_UCSC/hg19_UCSC.gtf
	annotationFile=/SAN/biomed/biomed14/vyp-scratch/temp_reference/hg19_UCSC/hg19_UCSC.gtf
	gffFile=/SAN/biomed/biomed14/vyp-scratch/temp_reference/hg19_UCSC/hg19_UCSC.gtf
	STARdir=/SAN/biomed/biomed14/vyp-scratch/temp_reference/hg19_UCSC/STAR
	;;
    human_hg38)
        fasta=${bigFilesBundleFolder}/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
        gtfFile=${bigFilesBundleFolder}/RNASeq/Human_hg38/Homo_sapiens.GRCh38.82_fixed.gtf
        gffFile=${bigFilesBundleFolder}/RNASeq/Human_hg38/Homo_sapiens.GRCh38.82_fixed.gff
        STARdir=${bigFilesBundleFolder}/RNASeq/Human_hg38/STAR 
	annotationFile=${bigFilesBundleFolder}/RNASeq/Human_hg38/biomart_annotations_human.tab
        #annotationFile=/scratch2/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/biomart_annotations_human.tab
	;;
    humanmuscle)
        refFolder=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Homo_sapiens/NCBI/build37.2
        IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome
        gtfFile=${refFolder}/Annotation/Genes/genes.gtf	
        gffFile=/cluster/project8/vyp/vincent/data/reference_genomes/gff/humanmuscle_iGenomes_NCBI37_with_ensembl.gff
        #annotationFile=${bigFilesBundleFolder}/human/biomart/biomart_annotations_human.tab
        annotationFile=/scratch2/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/biomart_annotations_human.tab
        geneModel=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/homoSapiens_geneTable_hg19_nochr.bed
        geneModelSummaryStats=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/homoSapiens_geneTable_hg19_chr1.bed
        db=hsapiens_gene_ensembl
        ;;
    Dict_Disc_masked)
        IndexBowtie2=${bigFilesBundleFolder}/RNASeq/Dict/dicty_masked_ERCC92
        gtfFile=${bigFilesBundleFolder}/RNASeq/Dict/dict_no_spike.gtf
        gffFile=MISSING
        annotationFile=not_done_yet
        ;;
    Dict_Disc)
        IndexBowtie2=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Dictyostelium_discoideum/sequence/Dictyostelium_discoideum.dictybase.01.23.dna.genome
        gtfFile=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Dictyostelium_discoideum/GTF/Dictyostelium_discoideum.dictybase.01.23.gtf
        gffFile=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Dictyostelium_discoideum/GTF/Dictyostelium_discoideum.dictybase.01.23.gff
        annotationFile=not_done_yet
        ;;
    pig)
        refFolder=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Sus_scrofa/NCBI/Sscrofa10.2
        IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome
        gtfFile=${refFolder}/Annotation/Genes/genes.gtf
        gffFile=${RNASEQBUNDLE}/pig/GTF/pig_iGenomes_NCBI_10_2_with_ensembl.gff
        annotationFile=${RNASEQBUNDLE}/pig/biomart/biomart_annotations_pig.tab
        ;;
    chicken)
        IndexBowtie2=${bigFilesBundleFolder}/RNASeq/Chicken/Gallus_gallus.Galgal4.dna.toplevel
        gtfFile=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/chicken/GTF/Gallus_gallus.Galgal4.78.gtf
        gffFile=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/chicken/GTF/Gallus_gallus.Galgal4.78.gff
        annotationFile=${RNASEQBUNDLE}/chicken/biomart/biomart_annotations_chicken.tab
        ;;
    rat)
        IndexBowtie2=${bigFilesBundleFolder}/RNASeq/Rat/Rattus_norvegicus.Rnor_5.0.dna_rm.toplevel
        gtfFile=${bigFilesBundleFolder}/RNASeq/Rat/Rattus_norvegicus.Rnor_5.0.79.gtf
        gffFile=${bigFilesBundleFolder}/RNASeq/Rat/Rattus_norvegicus.Rnor_5.0.79.gff
        annotationFile=${RNASEQBUNDLE}/rat/biomart/biomart_annotations_rat.tab
        ;;
    sheep)
        IndexBowtie2=${bigFilesBundleFolder}/RNASeq/Sheep/Ovis_aries.Oar_v3.1.dna_rm.toplevel
        gtfFile=${bigFilesBundleFolder}/RNASeq/Sheep/Ovis_aries.Oar_v3.1.80.gtf
        gffFile=${bigFilesBundleFolder}/RNASeq/Sheep/Ovis_aries.Oar_v3.1.80.gff
        annotationFile=${RNASEQBUNDLE}/sheep/biomart/biomart_annotations_sheep.tab
        ;;
    fly)
    	STARdir=${bigFilesBundleFolder}/RNASeq/Fly/STAR
	annotationFile=${bigFilesBundleFolder}/RNASeq/Fly/biomart_annotations_fly.tab
	gffFile=${bigFilesBundleFolder}/RNASeq/Fly/Drosophila_melanogaster.BDGP6.82.chr.corrected.names.gtf
	gtfFile=${bigFilesBundleFolder}/RNASeq/Fly/Drosophila_melanogaster.BDGP6.82.chr.corrected.names.gff
        #refFolder=${bigFilesBundleFolder}/RNASeq/Drosophila/Drosophila_melanogaster/NCBI/build5.41
        #IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome	
        #gtfFile=${RNASEQBUNDLE}/drosophila/GTF/Drosophila_melanogaster.BDGP5.75.gtf
        #gffFile=${RNASEQBUNDLE}/drosophila/GTF/Drosophila_melanogaster.BDGP5.75.gff
        #annotationFile=${RNASEQBUNDLE}/drosophila/biomart/biomart_annotations_drosophila.tab
        ;;
    dog)
        refFolder=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Canis_familiaris/NCBI/build3.1
        IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome	
        gtfFile=${refFolder}/Annotation/Genes/genes.gtf
        gffFile=${RNASEQBUNDLE}/dog/GTF/dog_iGenomes_NCBI_3_1_with_ensembl.gff
        annotationFile=${RNASEQBUNDLE}/dog/biomart/biomart_annotations_dog.tab
        ;;
    mouse)
        STARdir=${bigFilesBundleFolder}/RNASeq/Mouse/STAR
        annotationFile=${bigFilesBundleFolder}/RNASeq/Mouse/biomart_annotations_mouse.tab
        gffFile=${bigFilesBundleFolder}/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gff
        gtfFile=${bigFilesBundleFolder}/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gtf
        
	if [[ "$summary" == "refseq" ]];then
		gffFile=/SAN/vyplab/HuRNASeq/reference_datasets/mm10_refseq_genes_fixed.gff
		gtfFile=/SAN/vyplab/HuRNASeq/reference_datasets/mm10_refseq_genes_fixed.gtf
	fi
	;;
    tc1_mouse)
        refFolder=/SAN/biomed/biomed14/vyp-scratch/Zanda_AD_Tc1J20_RNASeq/Zanda_Tc1_reference/build1 
        IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome
        gtfFile=${RNASEQBUNDLE}/Tc1_mouse/GTF/Tc1.gtf #${refFolder}/Annotation/Genes/genes.gtf
        gffFile=${RNASEQBUNDLE}/Tc1_mouse/GTF/Tc1.gff  #${refFolder}/gff/tc1.gff
        cleanGtfFile=${RNASEQBUNDLE}/Tc1_mouse/GTF/Tc1.gtf #${refFolder}/Annotation/Genes/genes.gtf
        annotationFile=${RNASEQBUNDLE}/Tc1_mouse/tc1_annotations.tab #${refFolder}/annotations/biomart/tc1.tab
        geneModelSummaryStats=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/mm9_NCBI37_Ensembl_chr1.bed
        geneModel=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/mm9_NCBI37_Ensembl_nochr.bed
        ##db=mmusculus_gene_ensembl
        ;;
    *)
        stop "unknown species $species"
esac

files_exist $gtfFile $annotationFile $gffFile $STARdir

############### checking the input dataframe
h1=`awk '{if (NR == 1) print $1}'  $dataframe` 
h2=`awk '{if (NR == 1) print $2}'  $dataframe` 
h3=`awk '{if (NR == 1) print $3}'  $dataframe` 
if [[ "$h1" != "sample" ]]; then echo "header 1 must be sample"; exit; fi
if [[ "$h2" != "f1" ]]; then echo "header 2 must be f1 for fastq1"; exit; fi
if [[ "$h3" != "f2" ]]; then echo "header 3 must be f2 for fastq2"; exit; fi



hold=""
SCRATCH_DIR=/scratch0/${USER}/RNASeq_${code}
JAVA_DIR=${SCRATCH_DIR}/javastar 
# fastQC
function step0_QC {
   step0=${oFolder}/cluster/submission/step0.sh
	echo "
#$ -S /bin/bash
#$ -l h_vmem=4G,tmem=4G
#$ -l h_rt=72:00:00
#$ -pe smp 1
#$ -R y
#$ -o ${oFolder}/cluster/out
#$ -e ${oFolder}/cluster/error
#$ -N step0_${code}
#$ -wd ${oFolder}
echo \$HOSTNAME >&2
date >&2
" > $step0

	# create FASTQC output folder
	fastqcFolder=${oFolder}/fastqc/
	mkdir -p $fastqcFolder
	# read each row of the support table		
	tail -n +2  $dataframe | 
	while read sample f1 f2 condition;do
		# f1 
		f1array=(`echo $f1 | tr "," " "`)	# create bash array
        	f1_array_length=$(( ${#f1array[@]} - 1 )) # because bash arrays are 0 based

        	for i in `seq 0 ${f1_array_length}`;do
                	f1array[i]=${iFolder}/${f1array[i]} # append full path to each element of the array
                	files_exist ${f1array[i]} # check for existence
        		echo "
$fastqc -o ${fastqcFolder} ${f1array[i]} 
rename ${fastqcFolder}/${f1array[i]} ${fastqcFolder}/${sample}_${f1array[i]} ${fastqcFolder}/${f1array[i]}*
			" >> $step0
		done
		#f2
		if [[ "$f2" != "NA" ]];then
			f2array=(`echo $f2 | tr "," " "`)
			for i in `seq 0 ${f1_array_length}`;do
                		f2array[i]=${iFolder}/${f2array[i]} # append full path to each element of the array
                		files_exist ${f2array[i]} # check for existence
                	echo "
$fastqc -o $fastqcFolder ${f2array[i]} 
                	" >> $step0
        		done
		fi
	done
}

####################
# STEP1A: alignment
####################
function starSubmissionStep1a {
    starSubmissionStep1a=${oFolder}/cluster/submission/starSubmissionStep1a.sh
    STARoutput="BAM Unsorted"
    if [[ "$force" == "SJsOnly" ]];then
	STARoutput="None"
    fi
    echo "
#$ -S /bin/bash
#$ -l h_vmem=10G,tmem=10G
#$ -l h_rt=72:00:00
#$ -pe smp 4
#$ -R y
#$ -o ${oFolder}/cluster/out
#$ -e ${oFolder}/cluster/error
#$ -N step1a_${code}
#$ -l tscratch=10G  
#$ -wd ${oFolder}
echo \$HOSTNAME >&2
date >&2
mkdir -p $JAVA_DIR
" > $starSubmissionStep1a

# is trimming required?
if [ "$trim_galore" == "yes" ];then
	echo "
mkdir ${SCRATCH_DIR}/trimmed # make a scratch0 trimmed folder
" >> $starSubmissionStep1a
fi    

# for each sample in the support file
tail -n +2  $dataframe | 
while read sample f1 f2 condition;do
    if [[ "$f2" == "NA" ]]; then 
        paired=no 
    else paired=yes
    fi;
    echo "Sample $sample"
    finalOFolder=${oFolder}/${sample}
    dexseqfolder=${oFolder}/${sample}/dexseq
    mkdir -p ${finalOFolder} ${dexseqfolder}
    # go no further
    #if [[  -e ${finalOFolder}/${sample}_unique.bam.bai ]]
    #then
    #    echo ${finalOFolder}/${sample}_unique.bam.bai exists will not align again
    #    return 0
    #fi
   f1array=(`echo $f1 | tr "," " "`)
   f1_array_length=$(( ${#f1array[@]} - 1 )) # because bash arrays are 0 based
   for i in `seq 0 ${f1_array_length}`;do
   #	f1array[i]=${iFolder}/${f1array[i]} # append full path
	  files_exist ${iFolder}/${f1array[i]} # check for existence
   done

   echo Step1a: Aligning with STAR
   echo "paired = $paired"
# If sample is paired end        
    if [[ "$paired" == "yes" ]];then
        # create bash array of samples
        f2array=(`echo $f2 | tr "," " "`) 
        for i in `seq 0 ${f1_array_length}`;do
       #check for existence
	       files_exist ${iFolder}/${f2array[i]} 
        done
    
        # if Trim Galore is wanted run on paired end sample
        if [[ "$trim_galore" == "yes" ]];then
            echo $summary
    		if [[ "$summary" != "trimmed_exist" ]];then
    # trim each pair of files in the two arrays - assume that forward and reverse reads are in equal numbers of pieces
    			if [ ! -e ${iFolder}/trimmed ];then 
                    mkdir ${iFolder}/trimmed
                fi 
                # make trimmed folder
    			for i in `seq 0 $f1_array_length `;do 
                # i in length(array) - bash arrays are 0 based
    				echo "
# trim paired end with Trim Galore. Trim adapters and low quality (phred below 20). output to scratch0
$trimgalore --gzip -o ${SCRATCH_DIR}/trimmed --quality 20 --path_to_cutadapt $cutadapt --paired ${iFolder}/${f1array[i]} ${iFolder}/${f2array[i]} 
    "  >>  $starSubmissionStep1a
    			done
    		fi
    			#the trimmed files have a slightly different output
    			# prefix each fastq name with "trimmed/" so the trimmed outFolder is found
    		for i in `seq 0 $f1_array_length`;do
    			f1array[i]=`echo ${f1array[i]} | awk -F'/' '{print "trimmed/"$NF}' | sed 's/.fastq/_val_1.fq/g'`
    		done	
    		for i in `seq 0 $f1_array_length`;do
    			f2array[i]=`echo ${f2array[i]} | awk -F'/' '{print "trimmed/"$NF}' | sed 's/.fastq/_val_2.fq/g'`
    		done
    		echo ${f1array[@]}
    		echo ${f2array[@]}
    		if [[ "$summary" == "trimmed_exist" ]];then     
                            echo "
# trimmed_exist selected. Files have already been trimmed and exist in ${iFolder}/trimmed
" >> $starSubmissionStep1a      
                for i in `seq 0 $f1_array_length`;do
    				files_exist ${iFolder}/${f1array[i]}
    				files_exist ${iFolder}/${f2array[i]}
    			done
    		fi

                    #check that the trimming has happened. If not then exit
                    #echo "
    #if [ ! -e ${iFolder}/$f1 ]; then exit;fi" >> $starSubmissionStep1a
                #if QC step is wanted and ran successfully then the trimmed fastqs should be aligned.
    	fi

	# if the fastq files come in pieces, STAR takes them as A_R1.fastq,B_R1.fastq A_R2.fastq,B_R2.fastq. 
	# if files have just been trimmed, append scratch0. If files already trimmed or not trimmed then append iFolder. 
	# Separate by comma
	if [[ "$trim_galore" == "yes" && "$summary" != "trimmed_exist" ]];then
		fastqFolder=${SCRATCH_DIR}
	else
		fastqFolder=${iFolder}
	fi
	f1_total=`echo ${f1array[@]} | awk -v i=$fastqFolder 'BEGIN{RS=" ";ORS=","}{print i"/"$1}' | sed 's/,$//g'  `
	f2_total=`echo ${f2array[@]} | awk -v i=$fastqFolder 'BEGIN{RS=" ";ORS=","}{print i"/"$1}' | sed 's/,$//g'  `
	# two-pass mode - currently a hidden feature
    if [[ "$summary" == "twopass" ]]; then
        twopass="--twopassMode Basic"
    	memorymode="NoSharedMemory"
    else
        twopass=""
	memorymode="LoadAndKeep"
    fi

    # align with STAR	    
	echo "
# align with STAR. Output = ${STARoutput}
${starexec} --readFilesIn $f1_total $f2_total --readFilesCommand zcat --genomeLoad ${memorymode} --genomeDir ${STARdir} --runThreadN  4 --outSAMstrandField intronMotif --outFileNamePrefix ${SCRATCH_DIR}/${sample} --outSAMtype $STARoutput $twopass --outSAMunmapped Within --outSAMheaderHD ID:${sample} PL:Illumina
date >&2
# move the trimmed files back to trimmed folder in iFolder
mv -t ${iFolder}/trimmed `echo $f1_total | tr "," " " ` `echo $f2_total | tr "," 0" " `

# move the SJ.tab
mv ${SCRATCH_DIR}/${sample}SJ.out.tab ${finalOFolder}/

# move the unsorted bam file back to the outFolder for sorting in Step1b
mv ${SCRATCH_DIR}/${sample}Aligned.out.bam ${finalOFolder}/${sample}_unsorted.bam	

date >&2
# move all Log files
mv ${SCRATCH_DIR}/${sample}Log* ${finalOFolder}/
# remove old bam
rm ${SCRATCH_DIR}/${sample}Aligned.out.bam

	" >> $starSubmissionStep1a
#################### 
# IF SINGLE ENDED
####################
    elif [[ "$paired" == "no" ]]; then
	   fastqFolder=${iFolder}
       # if Trimming is desired, do it single ended
		
        if [[ "${trim_galore}" == "yes" ]];then
			echo $summary           
            fastqFolder=${SCRATCH_DIR}
			trimmedFolder=`dirname ${oFolder} | awk '{print $1"/trimmed/"}' `
			if [ ! -e $trimmedFolder ]; then mkdir ${trimmedFolder};fi
			
            if [[ "$summary" != "trimmed_exist" ]];then
                #trim each fastq separately
                for i in `seq 0 $f1_array_length `;do
					echo "
# trim single end with Trim Galore. Trim adapters and low quality (phred below 20)
$trimgalore --gzip -o ${SCRATCH_DIR}/trimmed --quality 20 --path_to_cutadapt $cutadapt ${iFolder}/${f1array[i]}
" >> $starSubmissionStep1a
				done
            fi
		#the trimmed files have a slightly different output
			for i in `seq 0 $f1_array_length`;do
				f1array[i]=`echo ${f1array[i]} | awk -F'/' '{print "trimmed/"$NF}' |  sed 's/.fastq/_trimmed.fq/g'`
        	done
        	if [[ "$summary" == "trimmed_exist" ]];then     
                fastqFolder=${iFolder}
				echo "
# trimmed_exist selected. Files have already been trimmed and exist in ${iFolder}/trimmed
" >> $starSubmissionStep1a      
         		for i in `seq 0 $f1_array_length`;do
         			files_exist ${iFolder}/${f1array[i]}
         		done
            fi
         fi

    # get fastq pieces in right format   
	f1_total=`echo ${f1array[@]} | 
    awk -v i=$fastqFolder 'BEGIN{RS=" ";ORS=","}{print i"/"$1}' | 
    sed 's/,$//g'  `
	# twopass mode
    if [[ "$summary" == "twopass" ]]; then
        twopass="--twopassMode Basic"
        memorymode="NoSharedMemory"
    else
        twopass=""
    memorymode="LoadAndKeep"
    fi

    # STAR alignment
    echo "
# align with STAR. Output = ${STARoutput}
${starexec} --readFilesIn ${f1_total} --readFilesCommand zcat --genomeLoad $memorymode $twopass --genomeDir ${STARdir} --runThreadN  4 --outSAMstrandField intronMotif --outFileNamePrefix ${SCRATCH_DIR}/${sample} --outSAMtype $STARoutput --outSAMunmapped Within --outSAMheaderHD ID:${sample} PL:Illumina
date >&2

# move the trimmed files back to trimmed folder in iFolder
mv -t ${iFolder}/trimmed `echo $f1_total | tr "," " " `

# move splice junction files
mv ${SCRATCH_DIR}/${sample}SJ.out.tab ${finalOFolder}/

# move the unsorted bam file back to the outFolder for sorting in Step1b
mv ${SCRATCH_DIR}/${sample}Aligned.out.bam ${finalOFolder}/${sample}_unsorted.bam  

date >&2

# move all Log files
mv ${SCRATCH_DIR}/${sample}Log* ${finalOFolder}/

# remove old bam file
#rm ${SCRATCH_DIR}/${sample}Aligned.out.bam

    " >> $starSubmissionStep1a
	fi 
done

echo "
# clean up 
rm -rf $JAVA_DIR
" >> $starSubmissionStep1a
echo "
# move trimming reports if created
if [ -e ${SCRATCH_DIR}/trimmed/*trimming* ];then
	mv ${SCRATCH_DIR}/trimmed/*trimming* ${iFolder}/trimmed/
fi

rm -rf ${SCRATCH_DIR}
" >> $starSubmissionStep1a

if [[ "$summary" == "twopass" ]];then
	echo "
# remove genome from memory
${starexec} --genomeLoad Remove --genomeDir ${STARdir}

" >> $starSubmissionStep1a
fi
}

########################################
# STEP1B: sorting and duplicate marking
########################################

function starSubmissionStep1b {
# per sample
# if the master table has been made before then remove it. Otherwise every time the submission script is run then more lines are appended to it.
  starMasterTableStep1b=${oFolder}/cluster/submission/starMasterTableStep1b.tab
  if [ -e $starMasterTableStep1b ]; then rm $starMasterTableStep1b;fi
  echo "scripts" > $starMasterTableStep1b
  tail -n +2  $dataframe | while read sample f1 f2 condition
  do
	finalOFolder=${oFolder}/${sample}
    SCRATCH_1b=/scratch0/${USER}/${sample}_sort/
    echo "
mkdir -p ${SCRATCH_1b}
# sort reads and mark duplicates with NovoSort. Write unique.bam back to original folder
$novosort --md --xs -f -t ${SCRATCH_1b} -6 -c 2 -m 12G ${finalOFolder}/${sample}_unsorted.bam -o ${finalOFolder}/${sample}_unique.bam

# if sorting is successful then remove the unsorted bam file
if [ -e ${finalOFolder}/${sample}_unique.bam ];then
    rm ${finalOFolder}/${sample}_unsorted.bam
fi
${samtools} index ${finalOFolder}/${sample}_unique.bam

${samtools} flagstat ${finalOFolder}/${sample}_unique.bam > ${finalOFolder}/${sample}_mappingStats.tab
" > ${oFolder}/cluster/submission/star_step1b_${sample}.sh
    echo "${oFolder}/cluster/submission/star_step1b_${sample}.sh" >> $starMasterTableStep1b
  done
  # the job array will start the SGE_TASK_ID at 2 and iterate through. 
  # Therefore the first 1b script will run at SGE_TASK_ID=2 
  # and the last will run at SGE_TASK_ID=(total number of scripts + 1)
  njobs1b=`wc -l $starMasterTableStep1b | awk '{print $1 + 1}'`
# submission script
  starSubmissionStep1b=${oFolder}/cluster/submission/starSubmissionStep1b.sh
  echo "
#$ -S /bin/bash
#$ -l h_vmem=7.75G
#$ -l tmem=7.75G
#$ -l h_rt=24:00:00
#$ -pe smp 2
#$ -R y
#$ -o ${oFolder}/cluster/out
#$ -e ${oFolder}/cluster/error
#$ -N step1b_${code}
#$ -l scr=20G
#$ -l tscratch=20G 
#$ -wd ${oFolder}
#$ -t 2-${njobs1b}
#$ -tc 20
echo \$HOSTNAME >&2
script=\`awk '{if (NR == '\$SGE_TASK_ID') print}' $starMasterTableStep1b\`
sh \$script
" > $starSubmissionStep1b
}

# dexseqCount
function starSubmissionStep2 {
# per sample
    starMasterTableStep2=${oFolder}/cluster/submission/starMasterTableStep2.tab
    echo "scripts" > $starMasterTableStep2

    #echo $PWD/${dataframe}; exit
    #echo ${oFolder}/${dataframe}; exit
    
    ##echo $dataframe; exit
    tail -n +2  $dataframe | while read sample f1 f2 condition; do
        if [[ "$f2" == "NA" ]]; then paired=no;  else paired=yes; fi;	
        dexseqfolder=${oFolder}/${sample}/dexseq

	finalOFolder=${oFolder}/${sample}
	echo "
$samtools view -F 0x0400 ${finalOFolder}/${sample}_unique.bam |  ${pythonbin} ${dexseqCount} --order=pos --paired=${paired} --stranded=${countStrand}  ${gffFile} - ${dexseqfolder}/${sample}_dexseq_counts.txt
$samtools view ${finalOFolder}/${sample}_unique.bam |  ${pythonbin} ${dexseqCount} --order=pos --paired=${paired} --stranded=${countStrand}  ${gffFile} - ${dexseqfolder}/${sample}_dexseq_counts_keep_dups.txt
" > ${oFolder}/cluster/submission/star_step2_${sample}.sh

        echo "${oFolder}/cluster/submission/star_step2_${sample}.sh" >> $starMasterTableStep2

    done
    #((nscripts=nscripts+1))
    njobs2=`wc -l $starMasterTableStep2 | awk '{print $1}'`
# submission script
    starSubmissionStep2=${oFolder}/cluster/submission/starSubmissionStep2.sh
    echo "
#$ -S /bin/bash
#$ -l h_vmem=5.9G
#$ -l tmem=5.9G
#$ -l h_rt=12:00:00
#$ -N step2_${code}
#$ -R y
#$ -o ${oFolder}/cluster/out
#$ -e ${oFolder}/cluster/error
#$ -wd ${oFolder}
#$ -t 2-${njobs2}
#$ -tc 20
echo \$HOSTNAME >&2
script=\`awk '{if (NR == '\$SGE_TASK_ID') print}' $starMasterTableStep2\`
sh \$script
" > $starSubmissionStep2
}




################################################# Now the scripts that take all samples together    

function starSubmissionStep3 {
# analyse all samples together
    starSubmissionStep3=${oFolder}/cluster/submission/starSubmissionStep3.sh    
    ncores=1
    nhours=0
    nminutes=0
    mem=0
    if [[ "$prepareCounts" == "yes" ]]; then mem=13.9; ((nminutes=nminutes+50)); fi
    if [[ "$Rdeseq" == "yes" ]]; then ((nhours=nhours+3)); mem=13.9; fi
#scale with the number of samples
    if [[ "$Rdexseq" == "yes" ]]; then ((nhours=nhours+18)); ncores=4;mem=5.9; fi
    #if [[ "$RpathwayGO" == "yes" ]]; then ((nhours=nhours+3)); mem=6; fi
    #if [[ "$RtopGO" == "yes" ]]; then ((nhours=nhours+3)); mem=6; fi
    echo "
#!/bin/bash
#$ -S /bin/bash
#$ -o ${clusterFolder}/out
#$ -e ${clusterFolder}/error
#$ -cwd
#$ -pe smp $ncores
#$ -l tmem=${mem}G,h_vmem=${mem}G
#$ -V
#$ -N step3_${code}
#$ -R y
#$ -l h_rt=${nhours}:${nminutes}:00
echo \$HOSTNAME >&2
" > $starSubmissionStep3
    if [[ "$prepareCounts" == "yes" ]]
    then
        files_exist ${annotationFile} ${dataframe} ${oFolder} ${countPrepareR}
    echo "
${Rscript} ${countPrepareR} --gff ${gffFile} --annotation.file ${annotationFile} --keep.dups ${keepDups} --support.frame ${dataframe} --code ${code} --iFolder ${oFolder} > ${clusterFolder}/R/count_prepare.out
    " >> $starSubmissionStep3
    fi
    ##############
    if [[ "$Rdeseq" == "yes" ]]
    then
        files_exist $deseqFinalProcessR $dataframe
    echo "
${Rscript} ${deseqFinalProcessR} --keep.sex ${keepSex} --support.frame ${dataframe} --keep.dups ${keepDups} --code ${code} --annotation.file ${annotationFile} --iFolder ${oFolder} > ${clusterFolder}/R/deseq_${stem}.out 
    " >> $starSubmissionStep3
    fi
    ##############
    if [[ "$Rdexseq" == "yes" ]]
    then
        files_exist $dexseqFinalProcessR $dataframe
    echo "
export LD_LIBRARY_PATH=/share/apps/zlib-1.2.8/lib:$LD_LIBRARY_PATH
${Rscript} ${dexseqFinalProcessR} --gff ${gffFile} --keep.sex ${keepSex} --keep.dups ${keepDups} --support.frame ${dataframe} --code ${code} --annotation.file ${annotationFile} --iFolder ${oFolder} > ${clusterFolder}/R/dexseq_${stem}.out
    " >> $starSubmissionStep3
    fi
    ##############
    if [[ "$RpathwayGO" == "yes" ]]
    then
        files_exist $pathwayGOAnalysisR $dataframe
    echo "
${Rscript} ${pathwayGOAnalysisR} --support.frame ${dataframe} --code ${code} --mart ${mart} --db ${db} --iFolder ${oFolder} > ${clusterFolder}/R/pathwayGO_${stem}.out 
    " >> $starSubmissionStep3
    fi
    ##############
    if [[ "$RtopGO" == "yes" ]]
    then
        files_exist $topGOAnalysisR $dataframe
    echo "
${Rscript} ${topGOAnalysisR} --support.frame ${dataframe} --code ${code} --mart ${mart} --db ${db} --iFolder ${oFolder} > ${clusterFolder}/R/topGO_${stem}.out 
    " >> $starSubmissionStep3
    fi
    ############# new GO analysis with GOSeq - very fast
    if [[ "$goseq" == "yes" ]];then
	files_exist $goseqR $dataframe
	echo "
${Rscript} ${goseqR} --species ${species} --oFolder ${oFolder} --mode DESeq
	" >> $starSubmissionStep3
	fi
}


### SUBMITTING TIME

if [[ "$step0_QC" == "yes" ]];then
    echo step0: fastQC
    step0_QC
    files_exist $step0
    if [[ "$submit" == "yes" ]];then
	qsub $step0
    fi	
fi



if [[ "$starStep1a" == "yes" ]]
then
   echo step1a: align
   starSubmissionStep1a
   files_exist $starSubmissionStep1a 
   if [[ "$submit" == "yes" ]]
   then
       qsub $hold $starSubmissionStep1a
       if [[ "$hold" == "" ]]; then hold="-hold_jid step1a_${code}"; else hold="$hold,step1b_${code}"; fi
   fi
fi

if [[ "$starStep1b" == "yes" ]]
then
   echo step1b: sorting and duplication removal
   starSubmissionStep1b
   files_exist $starSubmissionStep1b
   if [[ "$submit" == "yes" ]]
   then
       qsub $hold $starSubmissionStep1b
       if [[ "$hold" == "" ]]; then hold="-hold_jid step1b_${code}"; else hold="$hold,step1b_${code}"; fi
   fi
fi

if [[ "$starStep2" == "yes" ]]
then
   echo step2: dexseq count
   starSubmissionStep2
   files_exist $starSubmissionStep2
   if [[ "$submit" == "yes" ]]
   then
       qsub $hold $starSubmissionStep2
       if [[ "$hold" == "" ]]; then hold="-hold_jid step2_${code}"; else hold="$hold,step2_${code}"; fi
   fi
fi

if [[ "$prepareCounts" == "yes" || "$Rdeseq" == "yes" || "$Rdexseq" == "yes" || "$RpathwayGO" == "yes" || "$RtopGO" == "yes" || "$goseq" == "yes" ]]
then
    echo "Taking all samples together"
    starSubmissionStep3
    ls -ltrh $starSubmissionStep3
    if [[ "$submit" == "yes" ]]
    then
        qsub $hold $starSubmissionStep3
    fi
fi


