#!/bin/bash 

until [ -z "$1" ]
do
    # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
        --support)
            shift
            support=$1;;
        --pipelineBaseDir)
            shift
            pipelineBaseDir=$1;;
        --code) 
            shift 
            code=$1;; 
        --outputDir) 
            shift 
            outputDir=$1;;
        --clusterDir) 
            shift
            clusterDir=$1;;
        --step) 
            shift 
            step=$1;;
        --caseCondition) 
            shift 
            caseCondition=$1;;
        --species) 
            shift 
            species=$1;;
        -* )
            stop "Unrecognized option: $1"
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
    echo $1
done


step1a=${pipelineBaseDir}/sgseq_step1a.R
step1b=${pipelineBaseDir}/sgseq_step1b.R
step2=${pipelineBaseDir}/sgseq_step2.R

script_step1a=${clusterDir}/submission/sgseq_step1a.sh 
script_step1b=${clusterDir}/submission/sgseq_step1b.sh 
script_step2a=${clusterDir}/submission/sgseq_step2a.sh 
script_step2b=${clusterDir}/submission/sgseq_step2b.sh 
Rscript=/share/apps/R-3.3.2/bin/Rscript

# make directories
for dir in $outputDir ${clusterDir}/out ${clusterDir}/error ${clusterDir}/R; do
	if [ ! -e $dir ];then
		mkdir -p $dir
	fi
done

if [ $step = "step1a" ]; then  
echo "  
#$ -S /bin/bash
#$ -l h_vmem=12G,tmem=12G
#$ -l h_rt=72:00:00
#$ -pe smp 8  
#$ -R y
#$ -o ${clusterDir}/out
#$ -e ${clusterDir}/error
#$ -cwd 


$Rscript --vanilla ${step1a} --support.tab $support --code $code --output.dir $outputDir --species $species
" > $script_step1a

echo $script_step1a
fi

if [ $step = "step1b" ]; then  
echo "  
#$ -S /bin/bash
#$ -l h_vmem=12G,tmem=12G
#$ -l h_rt=72:00:00
#$ -pe smp 8  
#$ -R y
#$ -o ${clusterDir}/out
#$ -e ${clusterDir}/error
#$ -cwd 


$Rscript --vanilla ${step1b} --support.tab $support --code $code --output.dir $outputDir --case.condition $caseCondition  
" > $script_step1b

echo $script_step1b

fi

if [ $step = "step2a" ]; then  
echo "  
#$ -S /bin/bash
#$ -l h_vmem=8G,tmem=8G
#$ -l h_rt=72:00:00
#$ -pe smp 1  
#$ -R y
#$ -o ${clusterDir}/out
#$ -e ${clusterDir}/error
#$ -cwd 

$Rscript --vanilla ${step2} --step ${step} --support.tab $support --code $code --output.dir $outputDir 
" > $script_step2a

echo $script_step2a

fi

if [ $step = "step2b" ]; then  
echo "  
#$ -S /bin/bash
#$ -l h_vmem=8G,tmem=8G
#$ -l h_rt=72:00:00
#$ -pe smp 1  
#$ -R y
#$ -o ${clusterDir}/out
#$ -e ${clusterDir}/error
#$ -cwd 

$Rscript --vanilla ${step2} --step ${step} --support.tab $support --code $code --output.dir $outputDir 
" > $script_step2b

echo $script_step2b

fi
