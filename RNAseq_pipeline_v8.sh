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
tophatbin=${software}/tophat-2.0.13.Linux_x86_64/tophat2
bowtie2Folder=${software}/bowtie2-2.2.6
samtoolsFolder=${software}/samtools-0.1.19
samtools1=${software}/samtools-1.2/samtools

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
	--tophat)
	    shift
	    tophat=$1;;
	--miso)
	    shift
	    miso=$1;;
	--prepareCounts)
	    shift
	    prepareCounts=$1;;
	--dexseqcounts)
	    shift
	    dexseqcounts=$1;;
	--sampleQC)
	    shift
	    sampleQC=$1;;
	--runCufflinks)
	    shift
	    runCufflinks=$1;;
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



indexTranscriptome=${IndexBowtie2}_${species}_transcriptome
if [ ! -e ${indexTranscriptome}.rev.2.bt2 ]; then

    echo "No transcriptome bowtie index found, will now create it, which can take a while (about 30 minutes, but done once per species)"
    echo "$indexTranscriptome"
    echo "Script is create_transcriptome_index.sh, running it now"

    echo "
export PATH=${bowtie2Folder}:\$PATH

${tophatbin} --transcriptome-index=${indexTranscriptome} -G ${gtfFile} ${IndexBowtie2}
" > create_transcriptome_index.sh

    sh create_transcriptome_index.sh
    exit
fi



############### checking the input dataframe
h1=`awk '{if (NR == 1) print $1}'  $dataframe` 
h2=`awk '{if (NR == 1) print $2}'  $dataframe` 
h3=`awk '{if (NR == 1) print $3}'  $dataframe` 
if [[ "$h1" != "sample" ]]; then echo "header 1 must be sample"; exit; fi
if [[ "$h2" != "f1" ]]; then echo "header 2 must be f1 for fastq1"; exit; fi
if [[ "$h3" != "f2" ]]; then echo "header 3 must be f2 for fastq2"; exit; fi


if [[ "$tophat" == "yes" || "$miso" == "yes" || "$dexseqcounts" == "yes" || "$runCufflinks" == "yes" || "$sampleQC" == "yes" ]]; then  ##in this case we proceed per sample
    awk '{if ($1 != "sample") print}'  $dataframe | while read sample f1 f2 condition; do
	sleep 0.2

	echo "Sample $sample"
	if [[ "$f2" == "NA" ]]; then 
	    paired=no
	    viewOption=""
	else 
	    paired=yes
	    viewOption="-F 0x0400"
	fi

	
	
	#### create the folders
	finalOFolder=${oFolder}/${sample}
	cufflinksFolder=${finalOFolder}/cufflinks

	for folder in ${oFolder} ${finalOFolder}; do
	    if [ ! -e $folder ]; then mkdir $folder; fi
	done
	dexseqfolder=${finalOFolder}/dexseq
	
	#################################
	tophatlocal=$tophat
	if [[ "$tophatlocal" == "yes" && "$force" == "no" && -e ${finalOFolder}/${sample}_unique.bam.bai ]]; then tophatlocal=no; fi
	dexseqcountslocal=$dexseqcounts
	if [[ "$dexseqcountslocal" == "yes" && "$force" == "no" && -e ${dexseqfolder}/${sample}_dexseq_counts.txt ]]; then dexseqcountslocal=no; fi
	runCufflinkslocal=$runCufflinks
	if [[ "$runCufflinkslocal" == "yes" && "$force" == "no" && -e ${cufflinksFolder}/cufflinks_basic/genes.fpkm_tracking ]]; then runCufflinkslocal=no; fi
	

        ################################### the per sample jobs requirements for RAM and nb of cores, sorted by increasing requirements
	nhours=0
	memory=7
	ncores=1
	queue=blades
	
	req=1
	tmp_req=80
	if [[ "$miso" == "yes" ]]; then ((nhours=nhours+72)); memory=3.9; fi
	if [[ "$sampleQC" == "yes" ]]; then ((nhours=nhours+10)); memory=7; fi
	#if [[ "${dexseqcountslocal}" == "yes" ]]; then ((nhours=nhours+15)); memory=2; ncores=8; fi ## was 7
	if [[ "${dexseqcountslocal}" == "yes" ]]; then ((nhours=nhours+15)); memory=7; ncores=1; fi ## was 7
	if [[ "$tophatlocal" == "yes" ]]; then ((nhours=nhours+72)); ncores=4; memory=3.4; queue=novoalign; fi  ##4 days allowed ##was 1.5
	if [[ "$runCufflinks" == "yes" ]]; then ((nhours=nhours+10)); ncores=4; memory=3; queue=novoalign; fi
		
	
	######################## Now create the script itself
	script=${clusterFolder}/submission/X${sample}_${code}_RNASeq.sh
	echo "
#!/bin/bash
#$ -S /bin/bash
#$ -o ${clusterFolder}/out
#$ -e ${clusterFolder}/error
#$ -cwd
#$ -pe smp $ncores
#$ -l tmem=${memory}G,h_vmem=${memory}G
#$ -V
#$ -l scr=${tmp_req}G
#$ -R y
#$ -l h_rt=${nhours}:00:00

PYTHONPATH=/cluster/project8/vyp/vincent/libraries/python/share/apps/python-2.7.1/lib/python2.7/site-packages:/cluster/project8/vyp/vincent/libraries/python/lib/python:\$PYTHONPATH

export PYTHONPATH

" > $script
	
	
	if [[ "$tophatlocal" == "yes" ]]; then 
	    memory2=10  ##if we run tophat we know that 12G are available, so 11G to sort is OK
	    
            SCRATCH_DIR=/scratch0/${sample}
            DATA_DIR=${SCRATCH_DIR}/data
            RESULTS_DIR=${SCRATCH_DIR}/results
            TMP_DIR=${SCRATCH_DIR}/tmp
            JAVA_DIR=${SCRATCH_DIR}/java
	    
	    
	    echo "
export PATH=${bowtie2Folder}:\$PATH

mkdir $SCRATCH_DIR 
mkdir $DATA_DIR
mkdir $RESULTS_DIR
mkdir $TMP_DIR
mkdir $JAVA_DIR

" >> $script

	    
    	    ###################### first read of the pair
	    fullfile1=""
	    first1=TRUE
	    for lfile in `echo $f1 | sed -e 's/,/ /g'`; do ##for multiple files
		ifastq=${iFolder}/${lfile}
		ls -ltrh ${ifastq}
		if [ ! -e $ifastq ]; then echo "File $ifastq does not exist"; fi
		fastq1=`basename ${ifastq}` 

		echo "cp $ifastq ${DATA_DIR}/$fastq1" >> $script

		if [[ "$first1" == "TRUE" ]]; then
		    first1=FALSE
		    fullfile1=${DATA_DIR}/${fastq1}
		else
		    fullfile1="${fullfile1},${DATA_DIR}/${fastq1}"
		fi
	    done
	    	    
	    ###################### second read of the pair
	    fullfile2=""
	    if [[ "$f2" == "NA" ]]; then 
		ltype="--segment-length $segmentLength"
		fastqinput="${fastqFiles[ 1 ]}"
	    else 
		ltype="-r 220 --library-type $libstrand --segment-length $segmentLength"
		
		first2=TRUE
		for lfile in `echo $f2 | sed -e 's/,/ /g'`; do ##for multiple files        
   		    ifastq=${iFolder}/${lfile}
		    ls -ltrh ${ifastq}
		    if [ ! -e $ifastq ]; then echo "File $ifastq does not exist"; fi
		   
		    fastq2=`basename ${ifastq}`
		    echo "cp $ifastq $DATA_DIR/$fastq2" >> $script
		    
    
		    if [[ "$first2" == "TRUE" ]]; then
			first2=FALSE
			fullfile2=${DATA_DIR}/${fastq2}
		    else
		      fullfile2="${fullfile2},${DATA_DIR}/${fastq2}"
		    fi
		done
	    fi
	    
	    echo "
 
${tophatbin} --keep-fasta-order --transcriptome-index=${indexTranscriptome} --rg-id ${sample} --rg-sample ${sample} --rg-platform Illumina --tmp-dir $TMP_DIR --no-coverage-search -o ${RESULTS_DIR} -p $ncores ${ltype}  ${IndexBowtie2} $fullfile1 $fullfile2

cp -r ${RESULTS_DIR}/* ${finalOFolder} 

rm -rf ${RESULTS_DIR} ${DATA_DIR}

${samtools1} index ${finalOFolder}/accepted_hits.bam

$java -Xmx9g -jar ${picardDup} TMP_DIR=${JAVA_DIR} ASSUME_SORTED=true REMOVE_DUPLICATES=FALSE INPUT=${finalOFolder}/accepted_hits.bam OUTPUT=${finalOFolder}/${sample}_unique.bam METRICS_FILE=${finalOFolder}/metrics_${sample}_unique.tab

rm ${finalOFolder}/accepted_hits.bam ${finalOFolder}/accepted_hits.bam.bai

${samtools1} index ${finalOFolder}/${sample}_unique.bam

${samtools1} flagstat ${finalOFolder}/${sample}_unique.bam > ${finalOFolder}/${sample}_stats.txt

${samtools1} flagstat ${finalOFolder}/unmapped.bam > ${finalOFolder}/${sample}_unmapped_stats.txt

# Tidy up the scratch files 
rm -rf $SCRATCH_DIR 

" >> $script
		
	fi

#########
	if [[ "$miso" == "yes" ]]; then 
	    
	    MISOcode=`basename ${misoindex}`

	    misoMain=${finalOFolder}/miso
	    misoFolder=${finalOFolder}/miso/${MISOcode}

	    for folder in ${misoMain} ${misoFolder}; do
		if [ ! -e $folder ]; then mkdir $folder; fi 
	    done


	    echo "
rm -rf ${misoFolder}/*

${pythonbin} ${misoRunEvents} --compute-genes-psi ${misoindex} ${finalOFolder}/${sample}_unique.bam --output-dir ${misoFolder} --read-len 40 --paired-end 250 15

${pythonbin} ${runMiso} --summarize-samples ${misoFolder} ${misoFolder}
" >> $script
	fi

#######

	if [[ "${dexseqcountslocal}" == "yes" ]]; then
	    echo "Counting step"
	    memory2=6 ##We should have at least 7G or RAM available here so 6 is fine
	    

	    for folder in $dexseqfolder; do
		if [ ! -e $folder ]; then mkdir $folder; fi 
	    done

	    if [ ! -e ${gffFile} ]; then echo "Missing gff file ${gffFile}"; exit; fi
	    
	    
	    ###handling of the stranddeness below, a tad tricky really
	    countStrand=no
	    if [[ "$libstrand" == "fr-firststrand" ]]; then
		countStrand=yes
		countStrandReverse=reverse
	    fi
	    
	    if [[ "$libstrand" == "fr-secondstrand" ]]; then
		countStrand=reverse
		countStrandReverse=yes
	    fi
	    
	    echo "

${samtoolsFolder}/samtools view -F 0x0400 ${finalOFolder}/${sample}_unique.bam |  ${pythonbin} ${dexseqCount} --order=pos --paired=${paired} --stranded=${countStrand}  ${gffFile} - ${dexseqfolder}/${sample}_dexseq_counts.txt

${samtoolsFolder}/samtools view ${finalOFolder}/${sample}_unique.bam |  ${pythonbin} ${dexseqCount} --order=pos --paired=${paired} --stranded=${countStrand}  ${gffFile} - ${dexseqfolder}/${sample}_dexseq_counts_keep_dups.txt

" >> $script

	    if [[ "$stranded" == "yes" ]]; then
		echo "
$novosort -n -f -t /scratch0/ -0 -c $ncores -m ${memory2}G ${finalOFolder}/${sample}_unique.bam | ${samtoolsFolder}/samtools view -F 0x0400 - | ${pythonbin} ${dexseqCount} --paired=${paired} --stranded=${countStrandReverse} ${gffFile} - ${dexseqfolder}/${sample}_dexseq_counts_antisense.txt

$novosort -n -f -t /scratch0/ -0 -c $ncores -m ${memory2}G ${finalOFolder}/${sample}_unique.bam | ${samtoolsFolder}/samtools view - | ${pythonbin} ${dexseqCount} --paired=${paired} --stranded=${countStrandReverse} ${gffFile} - ${dexseqfolder}/${sample}_dexseq_counts_antisense_keep_dups.txt
" >> $script

	    fi
	fi


	if [[ "$sampleQC" == "yes" ]]; then

	QCfolder=${finalOFolder}/${sample}_QC
	if [ ! -e $QCfolder ]; then mkdir $QCfolder; fi
	
	    echo "

echo \"QC step\"

echo \"First where do reads map?\"
${pythonbin} ${rseqQCscripts}/read_distribution.py -i ${finalOFolder}/${sample}_unique.bam -r ${geneModelSummaryStats} >  ${QCfolder}/${sample}_summaryStats.tab 

echo \"Infer strandedness\"
${pythonbin} ${rseqQCscripts}/infer_experiment.py -i ${finalOFolder}/${sample}_unique.bam -r ${geneModel} >  ${QCfolder}/${sample}_inferExperiment.tab

echo \"Annotating junctions\"
${pythonbin} ${rseqQCscripts}/junction_annotation.py -i ${finalOFolder}/${sample}_unique.bam -r ${geneModel} -o ${QCfolder}/${sample}_splicing
${pythonbin} ${rseqQCscripts}/junction_saturation.py -i ${finalOFolder}/${sample}_unique.bam -r ${geneModel} -o ${QCfolder}/${sample}_splicing_saturation

echo \"Annotating read quality and duplication\"
${pythonbin} ${rseqQCscripts}/read_duplication.py -i ${finalOFolder}/accepted_hits.bam -o ${QCfolder}/${sample}_duplicate_prior_removal
${pythonbin} ${rseqQCscripts}/read_duplication.py -i ${finalOFolder}/${sample}_unique.bam -o ${QCfolder}/${sample}_duplicate_after_removal
${pythonbin} ${rseqQCscripts}/read_quality.py -i ${finalOFolder}/${sample}_unique.bam -o ${QCfolder}/${sample}_quality

echo \"Now computing RPKM within the QC step\"
${pythonbin} ${rseqQCscripts}/RPKM_count.py -r ${geneModel} -i ${finalOFolder}/${sample}_unique.bam -o ${QCfolder}/${sample}_RPKM

" >> $script

	fi

	
	if [[ "$runCufflinkslocal" == "yes" ]]; then   ##creates some RPKM data. Is it what we want to use? Not so sure.
	    echo "Cufflinks step"

	    cufflinksFolderBasic=${cufflinksFolder}/cufflinks_basic
	    cufflinksFolderBetter=${cufflinksFolder}/cufflinks_better

	    for folder in ${cufflinksFolder} ${cufflinksFolderBasic} ${cufflinksFolderBetter}; do
		if [ ! -e $folder ]; then mkdir $folder; fi
	    done

	    echo "

$cufflinks -p ${ncores} --library-type ${libstrand} -o ${cufflinksFolderBasic} --GTF $cleanGtfFile ${finalOFolder}/${sample}_unique.bam

#$cufflinks -p ${ncores} -o ${cufflinksFolderBetter} --GTF-guide $cleanGtfFile ${finalOFolder}/${sample}_unique.bam

" >> $script

	fi
	
	echo $script
	if [[ "$submit" == "yes" && "$nhours" != "0" ]]; then qsub $script; fi
	if [[ "$submit" == "local" ]]; then sh $script; fi
    done
fi




if [[ "$summary" == "yes" ]]; then
    
    outSum=${oFolder}/summary_reads.tab

    if [[ "$summaryRegions" == "yes" ]]; then
	extra=""
	for locregion in `echo $regionsCount | sed -e 's/,/ /g'`; do
	    extra="$extra\t$locregion"
	    echo -e "sample\tmappedReads\tdupReads\tunmappedReads\tRead1\tRead2$extra" > $outSum
	done
    else 
	echo -e "sample\tmappedReads\tdupReads\tunmappedReads\tRead1\tRead2" > $outSum
    fi

    

    awk '{if ($1 != "sample") print}' $dataframe | while read sample f1 f2 condition; do
	sFile=${oFolder}/${sample}/${sample}_stats.txt
	uFile=${oFolder}/${sample}/${sample}_unmapped_stats.txt
	echo "Looking at $sFile"
	ls -ltrh $sFile
	
	mappedReads=`awk '{if (NR == 1) {print $1}}' $sFile`
	dupReads=`awk '{if (NR == 4) {print $1}}' $sFile`
	unmappedReads=`awk '{if (NR == 1) {print $1}}' $uFile`
	Read1=`awk '{if (NR == 7) {print $1}}' $sFile`
	Read2=`awk '{if (NR == 8) {print $1}}' $sFile`

	
	if [[ "$summaryRegions" == "yes" ]]; then
	    extra=""
	    for locregion in `echo $regionsCount | sed -e 's/,/ /g'`; do
		echo "Counting reads for $locregion"
		extra="$extra\t`${samtoolsFolder}/samtools view -c ${oFolder}/${sample}/${sample}_unique.bam $locregion`"

	    done
	    echo -e "$sample\t$mappedReads\t$dupReads\t$unmappedReads\t$Read1\t$Read2$extra" >> $outSum
	else
	    echo -e "$sample\t$mappedReads\t$dupReads\t$unmappedReads\t$Read1\t$Read2" >> $outSum
	fi
	

		
    done

    ls -ltrh $outSum
fi


if [[ "$star" == "yes" ]]; then

    SCRATCH_DIR=/scratch0/${sample}
    JAVA_DIR=${SCRATCH_DIR}/java
    
    
    starScript=cluster/submission/star_RNASeq.sh

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

	finalOFolder=${oFolder}/${sample}
	if [ ! -e ${finalOFolder} ]; then mkdir ${finalOFolder}; fi
	
	if [[ "$force" == "yes" || ! -e ${finalOFolder}/${sample}_unique.bam.bai ]]; then
	    
	    
	    echo "
${starexec} --readFilesIn $f1 $f2 --readFilesCommand zcat --genomeLoad LoadAndKeep --genomeDir ${STARdir} --runThreadN  4 --outFileNamePrefix ${finalOFolder}/${sample} --outSAMtype BAM Unsorted

### now need to sort
$novosort -f -t /scratch0/ -0 -c 4 -m 20G ${finalOFolder}/${sample}Aligned.out.bam -o ${finalOFolder}/${sample}.bam

$java -Xmx9g -jar ${picardDup} TMP_DIR=${JAVA_DIR} ASSUME_SORTED=true REMOVE_DUPLICATES=FALSE INPUT=${finalOFolder}/${sample}.bam OUTPUT=${finalOFolder}/${sample}_unique.bam METRICS_FILE=${finalOFolder}/metrics_${sample}_unique.tab

${samtools1} index ${finalOFolder}/${sample}_unique.bam

rm ${finalOFolder}/${sample}.bam ${finalOFolder}/${sample}Aligned.out.bam 
" >> $starScript
	fi
    done

    echo "
rm -rf $JAVA_DIR
" >> $starScript

    ls -ltrh $starScript
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

    





