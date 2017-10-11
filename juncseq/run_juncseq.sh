#!/bin/bash 
Rscript=/share/apps/R-3.3.2/bin/Rscript

until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
        --gff)
            shift
            gff=$1;;
        --basedir)
            shift
            basedir=$1;;
        --support)
            shift
            support=$1;;
        --submit)
            shift
            submit=$1;;
        --baseoutputdir)
            shift
            baseoutputdir=$1;;
        --code)
            shift 
            code=$1;;
        --juncseqscript) 
            shift
            juncseqscript=$1;;
        -* )
            stop "Unrecognized option: $1"
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
    echo $1
done

clusterFolder=${baseoutputdir}/cluster/
submissiondir=${clusterFolder}/submission/
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

script=${submissiondir}/${code}_run_juncseq.sh
echo "#!/bin/bash
#$ -S /bin/bash
#$ -o ${clusterFolder}/out
#$ -e ${clusterFolder}/error
#$ -cwd
#$ -l tmem=11.8G,h_vmem=11.8G
#$ -V
#$ -l scr=5G
#$ -pe smp 6
#$ -R y
#$ -l h_rt=5:00:00

export LD_LIBRARY_PATH=/share/apps/zlib-1.2.8/lib:$LD_LIBRARY_PATH
$Rscript ${juncseqscript} --support.tab $support --gff $gff --basedir $basedir --outputdir $outputdir 
" > $script

if [ $submit == "yes" ]; then 
   qsub $script 
fi
 
