## version v7

#computer=vanHeel
computer=CS

if [[ "$computer" == "CS" ]]; then
    software=/cluster/project8/vyp/vincent/Software
    pythonbin=/share/apps/python-2.7.1/bin/python2.7
    ##Rbin=/cluster/project8/vyp/vincent/Software/R-3.1.2/bin/R
    Rbin=/cluster/project8/vyp/vincent/Software/R-3.2.2/bin/R

    misoRunEvents=/cluster/project8/vyp/vincent/Software/misopy-0.4.9/misopy/run_events_analysis.py
    runMiso=/cluster/project8/vyp/vincent/Software/misopy-0.4.9/misopy/run_miso.py

    javaTemp2="/scratch2/vyp-scratch2/vincent/java_temp"
    javaTemp="TMP_DIR=${javaTemp2}"
    java=/share/apps/jdk1.7.0_45/bin/java
    
    dexseqCount=/cluster/project8/vyp/vincent/libraries/R/installed/DEXSeq/python_scripts/dexseq_count.py    
    bigFilesBundleFolder=/scratch2/vyp-scratch2/reference_datasets/
    if [ ! -e $bigFilesBundleFolder ]; then bigFilesBundleFolder=/cluster/scratch3/vyp-scratch2/reference_datasets/
    fi
fi




if [[ "$computer" == "vanHeel" ]]; then
    software=/data_n2/vplagnol/Software
    pythonbin=/software/additional/epd-7.3.1/bin/python
    Rbin=/data_n2/vplagnol/Software/R-3.0.2/bin/R
    
    dexseqCount=/data_n2/vplagnol/Rlibs/installed/DEXSeq/python_scripts/dexseq_count.py

    javaTemp2="/data_n1/vanheel_singlecellgenomics/tmp"
    javaTemp="TMP_DIR=${javaTemp2}"

fi


echo "Base of RNA-Seq pipeline is located here: $RNASEQPIPBASE"



countPrepareR=${RNASEQPIPBASE}/counts_prepare_pipeline.R
dexseqFinalProcessR=${RNASEQPIPBASE}/dexseq_pipeline_v2.R
deseqFinalProcessR=${RNASEQPIPBASE}/deseq2_pipeline.R
pathwayGOAnalysisR=${RNASEQPIPBASE}/pathwayGO_pipeline.R
topGOAnalysisR=${RNASEQPIPBASE}/topGO_pipeline.R
novosort=${software}/novocraft/novosort


for file in $countPrepareR; do
    if [ ! -e $file ]; then echo "Missing script fule $countPrepareR"; fi
done



starexec=/cluster/project8/vyp/vincent/Software/STAR/bin/Linux_x86_64_static/STAR
samtools=${software}/samtools-1.2/samtools

cufflinks=${software}/cufflinks-2.1.1.Linux_x86_64/cufflinks

rseqQCscripts=${software}/RSeQC-2.3.3/scripts

picardDup=${software}/picard-tools-1.100/MarkDuplicates.jar
picardStats=${software}/picard-tools-1.100/BamIndexStats.jar
picardMetrics=${software}/picard-tools-1.100/CalculateHsMetrics.jar
picardReorder=${software}/picard-tools-1.100/ReorderSam.jar


keepSex=FALSE  ## should sex chromosomes be kept in the differential expression analysis?
superLong=no
force=yes
species=mouse
segmentLength=25

mart=ensembl
db=mmusculus_gene_ensembl

summary=no
summaryRegions=no
prepareCounts=no
Rdexseq=no
Rdeseq=no
RpathwayGO=no
RtopGO=no
oFolder=temp ##default output


for folder in $oFolder; do
    if [ ! -e $folder ]; then mkdir $folder; fi
done


tophat=no
miso=no
runCufflinks=no
submit=yes

stranded=no
libstrand=fr-unstranded
keepDups=FALSE

code=""
dataframe=none
iFolder=""
misoindex=NA
oFolder=RNAseq_processed
mainscript=combined.sh
stem=""

until [ -z "$1" ]; do
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
	--star)
	    shift
	    star=$1;;
	--superLong)
	    shift
	    superLong=$1;;
	--keep.sex)
	    shift
	    keepSex=$1;;
	--misoindex)
	    shift
	    misoindex=$1;;
	--summary)
	    shift
	    summary=$1;;
	--summaryRegions)
	    shift
	    summaryRegions=yes
	    regionsCount=$1;;
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
	--mainscript)
            shift
            mainscript=$1;;
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
	--segmentLength)
            shift
            segmentLength=$1;;
	--submit)
	    shift
            submit=$1;;
        --stem)
            shift
            stem=$1;;
	--keepDups)
	    keepDups=TRUE
	    shift;;
	-* )
	    echo "Unrecognized option: $1"
	    exit 1;;
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
    echo $1
done


echo "Strand information $stranded $libstrand"
########################## estimate the duration of the jobs


if [[ "$stem" == "" ]]; then stem=$code; fi

if [[ "$superLong" == "yes" ]]; then ((nhours=nhours+nhours)); fi


## create the output folders
clusterFolder=${oFolder}/cluster

for folder in ${oFolder} ${clusterFolder} ${clusterFolder}/out ${clusterFolder}/error ${clusterFolder}/R  ${clusterFolder}/submission; do
    if [ ! -e $folder ]; then mkdir $folder; fi
done


if [ ! -e $dataframe ]; then
    echo "File $dataframe is required but does not exist"
fi
cp "$dataframe" "${oFolder}/"


if [[ "$species" == "zebrafish" ]]; then
    refFolder=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Danio_rerio/NCBI/Zv9
    gtfFile=${refFolder}/Annotation/Genes/genes.gtf    
    fasta=${refFolder}/Sequence/WholeGenomeFasta/genome.fa
    IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome

    gffFile=${RNASEQBUNDLE}/zebrafish/GTF/zebrafish_iGenomes_Zv9_with_ensembl.gff
    annotationFile=${RNASEQBUNDLE}/zebrafish/biomart/biomart_annotations_zebrafish.tab

    #geneModel=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/zebraFish_refSeqTable_zv9_nochr.bed
    #geneModelSummaryStats=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/zebraFish_refSeqTable_zv9_chr1.bed
fi

if [[ "$species" == "DvH_sc_human" ]]; then
    refFolder=/data_n1/vanheel_singlecellgenomics/support/Homo_sapiens/NCBI/build37.2
    gtfFile=${refFolder}/Annotation/Genes/genes.gtf    
    fasta=${refFolder}/Sequence/WholeGenomeFasta/genome.fa
    IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/human_b37_with_spikes
fi


if [[ "$species" == "human" ]]; then
    refFolder=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Homo_sapiens/NCBI/build37.2
    gtfFile=${refFolder}/Annotation/Genes/genes.gtf    
    fasta=${refFolder}/Sequence/WholeGenomeFasta/genome.fa
    IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome

    ### stuff below should go to the bundle
    SNPlist=${RNASEQBUNDLE}/human/exonic_snps_human_hg19_clean.tab
    gffFile=${RNASEQBUNDLE}/human/GTF/human_iGenomes_NCBI37_with_ensembl.gff
    cleanGtfFile=${RNASEQBUNDLE}/human/GTF/human_iGenomes_NCBI37_with_ensembl.gtf
    annotationFile=${RNASEQBUNDLE}/human/biomart/biomart_annotations_human.tab
    geneModel=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/homoSapiens_geneTable_hg19_nochr.bed
    geneModelSummaryStats=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/homoSapiens_geneTable_hg19_chr1.bed

    db=hsapiens_gene_ensembl
fi


if [[ "$species" == "human_hg38" ]]; then
    if [ -e /cluster/scratch3/vyp-scratch2/ ]; then
	    refFolder=/cluster/scratch3/vyp-scratch2/reference_datasets    
	else 
	    refFolder=/scratch2/vyp-scratch2/reference_datasets   
	fi
    fasta=${bigFilesBundleFolder}/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
    IndexBowtie2=${bigFilesBundleFolder}/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
    
    gffFile=${bigFilesBundleFolder}/human_reference_sequence/GTF_files/Homo_sapiens_GRCh38_78_fixed.gff
    gtfFile=${bigFilesBundleFolder}/human_reference_sequence/GTF_files/Homo_sapiens_GRCh38_78_fixed.gtf

    annotationFile=${RNASEQBUNDLE}/human_hg38/biomart/biomart_annotations_human.tab

fi

if [[ "$species" == "humanmuscle" ]]; then
    refFolder=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Homo_sapiens/NCBI/build37.2
    IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome
    gtfFile=${refFolder}/Annotation/Genes/genes.gtf	

    gffFile=/cluster/project8/vyp/vincent/data/reference_genomes/gff/humanmuscle_iGenomes_NCBI37_with_ensembl.gff
    annotationFile=${RNASEQBUNDLE}/human/biomart/biomart_annotations_human.tab
    geneModel=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/homoSapiens_geneTable_hg19_nochr.bed
    geneModelSummaryStats=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/homoSapiens_geneTable_hg19_chr1.bed

    db=hsapiens_gene_ensembl
fi



if [[ "$species" == "Dict_Disc_masked" ]]; then
    
    IndexBowtie2=${bigFilesBundleFolder}/RNASeq/Dict/dicty_masked_ERCC92
    gtfFile=${bigFilesBundleFolder}/RNASeq/Dict/dict_no_spike.gtf
    gffFile=MISSING

    annotationFile=not_done_yet

fi



if [[ "$species" == "Dict_Disc" ]]; then
    
    IndexBowtie2=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Dictyostelium_discoideum/sequence/Dictyostelium_discoideum.dictybase.01.23.dna.genome
    gtfFile=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Dictyostelium_discoideum/GTF/Dictyostelium_discoideum.dictybase.01.23.gtf
    gffFile=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Dictyostelium_discoideum/GTF/Dictyostelium_discoideum.dictybase.01.23.gff

    annotationFile=not_done_yet

fi

if [[ "$species" == "pig" ]]; then

    refFolder=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Sus_scrofa/NCBI/Sscrofa10.2
    IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome
    gtfFile=${refFolder}/Annotation/Genes/genes.gtf

    gffFile=${RNASEQBUNDLE}/pig/GTF/pig_iGenomes_NCBI_10_2_with_ensembl.gff
    annotationFile=${RNASEQBUNDLE}/pig/biomart/biomart_annotations_pig.tab
fi

if [[ "$species" == "chicken" ]]; then
    IndexBowtie2=${bigFilesBundleFolder}/RNASeq/Chicken/Gallus_gallus.Galgal4.dna.toplevel
    gtfFile=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/chicken/GTF/Gallus_gallus.Galgal4.78.gtf
    gffFile=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/chicken/GTF/Gallus_gallus.Galgal4.78.gff
    annotationFile=${RNASEQBUNDLE}/chicken/biomart/biomart_annotations_chicken.tab
fi


if [[ "$species" == "rat" ]]; then
    IndexBowtie2=${bigFilesBundleFolder}/RNASeq/Rat/Rattus_norvegicus.Rnor_5.0.dna_rm.toplevel
    gtfFile=${bigFilesBundleFolder}/RNASeq/Rat/Rattus_norvegicus.Rnor_5.0.79.gtf
    gffFile=${bigFilesBundleFolder}/RNASeq/Rat/Rattus_norvegicus.Rnor_5.0.79.gff
    annotationFile=${RNASEQBUNDLE}/rat/biomart/biomart_annotations_rat.tab
fi


if [[ "$species" == "sheep" ]]; then
    IndexBowtie2=${bigFilesBundleFolder}/RNASeq/Sheep/Ovis_aries.Oar_v3.1.dna_rm.toplevel
    gtfFile=${bigFilesBundleFolder}/RNASeq/Sheep/Ovis_aries.Oar_v3.1.80.gtf
    gffFile=${bigFilesBundleFolder}/RNASeq/Sheep/Ovis_aries.Oar_v3.1.80.gff
    annotationFile=${RNASEQBUNDLE}/sheep/biomart/biomart_annotations_sheep.tab
fi


if [[ "$species" == "drosophila" ]]; then
    refFolder=${bigFilesBundleFolder}/RNASeq/Drosophila/Drosophila_melanogaster/NCBI/build5.41
    IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome	
    gtfFile=${RNASEQBUNDLE}/drosophila/GTF/Drosophila_melanogaster.BDGP5.75.gtf
    gffFile=${RNASEQBUNDLE}/drosophila/GTF/Drosophila_melanogaster.BDGP5.75.gff
    annotationFile=${RNASEQBUNDLE}/drosophila/biomart/biomart_annotations_drosophila.tab
fi


if [[ "$species" == "dog" ]]; then
    refFolder=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Canis_familiaris/NCBI/build3.1
    IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome	
    gtfFile=${refFolder}/Annotation/Genes/genes.gtf

    gffFile=${RNASEQBUNDLE}/dog/GTF/dog_iGenomes_NCBI_3_1_with_ensembl.gff
    annotationFile=${RNASEQBUNDLE}/dog/biomart/biomart_annotations_dog.tab

fi

if [[ "$species" == "mouse" ]]; then
    STARdir=${bigFilesBundleFolder}/RNASeq/Mouse/STAR
    annotationFile=${bigFilesBundleFolder}/RNASeq/Mouse/biomart_annotations_mouse.tab
    gffFile=${bigFilesBundleFolder}/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gff
    gtfFile=${bigFilesBundleFolder}/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gtf
fi



if [[ "$species" == "tc1_mouse" ]]; then

    refFolder=/SAN/biomed/biomed14/vyp-scratch/Zanda_AD_Tc1J20_RNASeq/Zanda_Tc1_reference/build1 
    IndexBowtie2=${refFolder}/Sequence/Bowtie2Index/genome
	
    gtfFile=${RNASEQBUNDLE}/Tc1_mouse/GTF/Tc1.gtf #${refFolder}/Annotation/Genes/genes.gtf
    gffFile=${RNASEQBUNDLE}/Tc1_mouse/GTF/Tc1.gff  #${refFolder}/gff/tc1.gff
    cleanGtfFile=${RNASEQBUNDLE}/Tc1_mouse/GTF/Tc1.gtf #${refFolder}/Annotation/Genes/genes.gtf
    annotationFile=${RNASEQBUNDLE}/Tc1_mouse/tc1_annotations.tab #${refFolder}/annotations/biomart/tc1.tab


    geneModelSummaryStats=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/mm9_NCBI37_Ensembl_chr1.bed
    geneModel=/cluster/project8/vyp/vincent/data/reference_genomes/gene_tables/mm9_NCBI37_Ensembl_nochr.bed
    
    ##db=mmusculus_gene_ensembl

fi



############### checking the input dataframe
h1=`awk '{if (NR == 1) print $1}'  $dataframe` 
h2=`awk '{if (NR == 1) print $2}'  $dataframe` 
h3=`awk '{if (NR == 1) print $3}'  $dataframe` 
if [[ "$h1" != "sample" ]]; then echo "header 1 must be sample"; exit; fi
if [[ "$h2" != "f1" ]]; then echo "header 2 must be f1 for fastq1"; exit; fi
if [[ "$h3" != "f2" ]]; then echo "header 3 must be f2 for fastq2"; exit; fi




if [[ "$star" == "yes" ]]; then

    SCRATCH_DIR=/scratch0/${sample}
    JAVA_DIR=${SCRATCH_DIR}/java
    
    
    starScript=cluster/submission/star_RNASeq.sh
    starMasterTableStep2=${oFolder}/cluster/submission/starMasterTableStep2.tab
    starMasterTableStep3=${oFolder}/cluster/submission/starMasterTableStep3.tab

    starSubmissionStep2=${oFolder}/cluster/submission/starSubmissionStep2.sh
    starSubmissionStep3=${oFolder}/cluster/submission/starSubmissionStep3.sh

    echo "scripts" > $starMasterTableStep2
    echo "scripts" > $starMasterTableStep3
    
    echo "#$ -S /bin/bash
#$ -l h_vmem=8.4G
#$ -l tmem=8.4G
#$ -l h_rt=12:00:00
#$ -pe smp 4
#$ -R y
#$ -o cluster/out
#$ -e cluster/error
#$ -cwd

mkdir $JAVA_DIR
" > $starScript




    tail -n +2  $dataframe | while read sample f1 f2 condition; do

	echo "Sample $sample"
	finalOFolder=${oFolder}/${sample}
	dexseqfolder=${oFolder}/${sample}/dexseq

	if [ ! -e ${finalOFolder} ]; then mkdir ${finalOFolder}; fi
	if [ ! -e ${dexseqfolder} ]; then mkdir ${dexseqfolder}; fi
	
	if [[ "$force" == "yes" || ! -e ${finalOFolder}/${sample}_unique.bam.bai ]]; then
	    
	    echo "
${starexec} --readFilesIn $f1 $f2 --readFilesCommand zcat --genomeLoad LoadAndKeep --genomeDir ${STARdir} --runThreadN  4 --outFileNamePrefix ${finalOFolder}/${sample} --outSAMtype BAM Unsorted
" > $starScript
	    
	    echo "
$novosort -f -t /scratch0/ -0 -c 1 -m 7G ${finalOFolder}/${sample}Aligned.out.bam -o ${finalOFolder}/${sample}.bam

$java -Xmx9g -jar ${picardDup} TMP_DIR=${JAVA_DIR} ASSUME_SORTED=true REMOVE_DUPLICATES=FALSE INPUT=${finalOFolder}/${sample}.bam OUTPUT=${finalOFolder}/${sample}_unique.bam METRICS_FILE=${finalOFolder}/metrics_${sample}_unique.tab

${samtools} index ${finalOFolder}/${sample}_unique.bam

rm ${finalOFolder}/${sample}.bam ${finalOFolder}/${sample}Aligned.out.bam 

" > ${oFolder}/cluster/submission/star_step2_${sample}.sh
	    

	    echo "
$samtools view -F 0x0400 ${finalOFolder}/${sample}_unique.bam |  ${pythonbin} ${dexseqCount} --order=pos --paired=${paired} --stranded=${countStrand}  ${gffFile} - ${dexseqfolder}/${sample}_dexseq_counts.txt

$samtools view ${finalOFolder}/${sample}_unique.bam |  ${pythonbin} ${dexseqCount} --order=pos --paired=${paired} --stranded=${countStrand}  ${gffFile} - ${dexseqfolder}/${sample}_dexseq_counts_keep_dups.txt

" > ${oFolder}/cluster/submission/star_step3_${sample}.sh


	    echo "${oFolder}/cluster/submission/star_step2_${sample}.sh" >> $starMasterTableStep2
	    echo "${oFolder}/cluster/submission/star_step3_${sample}.sh" >> $starMasterTableStep3
	    ((nscripts=nscripts+1))

	fi

    done

    njobs3=`wc -l $starMasterTableStep3 | awk '{print $1}'`

    echo "#$ -S /bin/bash
#$ -l h_vmem=8.4G
#$ -l tmem=8.4G
#$ -l h_rt=12:00:00
#$ -pe smp 1
#$ -R y
#$ -o ${oFolder}/cluster/out
#$ -e ${oFolder}/cluster/error
#$ -wd ${oFolder}
#$ -t 1-${nscripts}
#$ -tc 20

script=\`awk '{if (NR == '\$SGE_TASK_ID') print}' $starMasterTableStep2\`

sh \$script

" > $starSubmissionStep2

    echo "#$ -S /bin/bash
#$ -l h_vmem=8.4G
#$ -l tmem=8.4G
#$ -l h_rt=12:00:00
#$ -pe smp 1
#$ -R y
#$ -o ${oFolder}/cluster/out
#$ -e ${oFolder}/cluster/error
#$ -wd ${oFolder}
#$ -t 2-${njobs3}
#$ -tc 20

script=\`awk '{if (NR == '\$SGE_TASK_ID') print}' $starMasterTableStep3\`

sh \$script

" > $starSubmissionStep3


    echo "
rm -rf $JAVA_DIR
" >> $starScript

    ls -ltrh $starScript $starSubmissionStep2 $starSubmissionStep3
fi


################################################# Now the scripts that take all samples together    

if [[ "$prepareCounts" == "yes" || "$Rdeseq" == "yes" || "$Rdexseq" == "yes" || "$RpathwayGO" == "yes" || "$RtopGO" == "yes" ]]; then
    

    

    ncores=1
    nhours=0
    nminutes=0
    mem=0
    
    
    if [[ "$prepareCounts" == "yes" ]]; then mem=13.9; ((nminutes=nminutes+50)); fi
    if [[ "$Rdeseq" == "yes" ]]; then ((nhours=nhours+3)); mem=13.9; fi
    if [[ "$Rdexseq" == "yes" ]]; then ((nhours=nhours+18)); ncores=4;mem=5.9; fi
    if [[ "$RpathwayGO" == "yes" ]]; then ((nhours=nhours+3)); mem=6; fi
    if [[ "$RtopGO" == "yes" ]]; then ((nhours=nhours+3)); mem=6; fi

    echo "
#!/bin/bash
#$ -S /bin/bash
#$ -o ${clusterFolder}/out
#$ -e ${clusterFolder}/error
#$ -cwd
#$ -pe smp $ncores
#$ -l tmem=${mem}G,h_vmem=${mem}G
#$ -V
#$ -R y
#$ -l h_rt=${nhours}:${nminutes}:00

" > $mainscript
	

    if [[ "$prepareCounts" == "yes" ]]; then

	echo "

${Rbin} CMD BATCH --no-save --no-restore --gff=${gffFile} --annotation.file=${annotationFile} --keep.dups=${keepDups} --support.frame=${dataframe} --code=${code} --iFolder=${oFolder} ${countPrepareR} ${clusterFolder}/R/count_prepare.out

" >> $mainscript
	
    fi


##############
    if [[ "$Rdeseq" == "yes" ]]; then

        for file in $deseqFinalProcessR $dataframe; do
            if [ ! -e $file ]; then echo "$file does not exist"; exit; fi
	done
	

	echo "
${Rbin} CMD BATCH --no-save --no-restore --keep.sex=${keepSex} --support.frame=${dataframe} --keep.dups=${keepDups} --code=${code} --annotation.file=${annotationFile} --iFolder=${oFolder} ${deseqFinalProcessR} ${clusterFolder}/R/deseq_${stem}.out 
" >> $mainscript

    fi

##############
    if [[ "$Rdexseq" == "yes" ]]; then
	
	for file in $dexseqFinalProcessR $dataframe; do
	    if [ ! -e $file ]; then echo "$file does not exist"; exit; fi
	done

	echo "
${Rbin} CMD BATCH --no-save --no-restore --gff=${gffFile} --keep.sex=${keepSex} --keep.dups=${keepDups} --support.frame=${dataframe} --code=${code} --annotation.file=${annotationFile} --iFolder=${oFolder} ${dexseqFinalProcessR} ${clusterFolder}/R/dexseq_${stem}.out

" >> $mainscript
    fi


##############
    if [[ "$RpathwayGO" == "yes" ]]; then

        for file in $pathwayGOAnalysisR $dataframe; do
            if [ ! -e $file ]; then echo "$file does not exist"; exit; fi
	done

	echo "
${Rbin} CMD BATCH --no-save --no-restore --support.frame=${dataframe} --code=${code} --mart=${mart} --db=${db} --iFolder=${oFolder} ${pathwayGOAnalysisR} ${clusterFolder}/R/pathwayGO_${stem}.out 
" >> $mainscript

    fi
    

##############
    if [[ "$RtopGO" == "yes" ]]; then

        for file in $topGOAnalysisR $dataframe; do
            if [ ! -e $file ]; then echo "$file does not exist"; exit; fi
	done
	
	echo "
${Rbin} CMD BATCH --no-save --no-restore --support.frame=${dataframe} --code=${code} --mart=${mart} --db=${db} --iFolder=${oFolder} ${topGOAnalysisR} ${clusterFolder}/R/topGO_${stem}.out 
" >> $mainscript

    fi
    

#############
    
    echo $mainscript
	

fi

    




