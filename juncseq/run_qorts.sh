#!/bin/bash 

## Software 
JAVA=/usr/java/latest/bin/java
QoRTsJar=/cluster/project8/vyp/kitty/Software/QoRTs.jar

until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
        --gtf)
            shift
            gtf=$1;;
        --basedir)
            shift
            basedir=$1;;
        --support)
            shift
            support=$1;;
        --submit)
            shift
            submit=$1;;
        --clusterFolder)
            shift
            clusterFolder=$1;;
        --baseoutputdir)
            shift
            baseoutputdir=$1;;
        -* )
            stop "Unrecognized option: $1"
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
    echo $1
done

clusterFolder=${baseoutputdir}/cluster/
submissiondir=${clusterFolder}/submission 
clusteroutputdir=${clusterFolder}/out
errordir=${clusterFolder}/error 

outputdir=${baseoutputdir}/qorts/ 

## Check the correct directories exist, if not, make them 
if [ ! -e $clusterFolder ]; then 
   mkdir $clusterFolder 
fi 

if [ ! -e $submissiondir ]; then 
   mkdir $submissiondir 
fi 

if [ ! -e $clusteroutputdir ]; then 
   mkdir $clusteroutputdir 
fi 

if [ ! -e $errordir ]; then 
   mkdir $errordir 
fi

awk '{if ($1 != "sample") print}'  $support | while read sample f1 f2 condition; do
   qortsdir=${outputdir}/${sample}/
   if [ ! -e $qortsdir ]; then 
       mkdir -p $qortsdir 
   fi 
   script=${submissiondir}/${code}_${sample}_qorts.sh
   bam=${basedir}/${sample}/${sample}_unique.bam 
   echo $bam 
   if [ ! -e $errordir ]; then 
       echo "Bam file not found?!" 
   fi
echo "#!/bin/bash
#$ -S /bin/bash
#$ -o ${clusterFolder}/out
#$ -e ${clusterFolder}/error
#$ -cwd
#$ -l tmem=11.8G,h_vmem=11.8G
#$ -V
#$ -l scr=5G
#$ -R y
#$ -l h_rt=5:00:00

export MALLOC_ARENA_MAX=4
$JAVA -Xmx10G -XX:ParallelGCThreads=1 -jar $QoRTsJar  QC --stranded --runFunctions writeKnownSplices,writeNovelSplices,writeSpliceExon $bam $gtf $qortsdir 

" > $script

if [ $submit == "yes" ]; then 
   qsub $script 
fi
 
done

