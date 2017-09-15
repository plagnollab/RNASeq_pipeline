pipelineBaseDir=/SAN/vyplab/HuRNASeq/RNASeq_pipeline/SGSeq/
sgseq_script=${pipelineBaseDir}/run_sgseq.sh

support=/SAN/vyplab/HuRNASeq/Fly_C9/fly_SGSeq_support.tab
code=Fly_C9
species=fly
step=step1b # step1a - just annotated events; step1b - novel events too
outputDir=/SAN/vyplab/HuRNASeq/Fly_C9/SGSeq
submit=no
#support must be sample_name    file_bam        condition

if [ ! -e $outputDir ]; then
                mkdir $outputDir
fi


sh $sgseq_script --species ${species} \
                 --step ${step} \
                 --support ${support} \
                 --code ${code} \
                 --submit ${submit} \
                 --outputDir ${outputDir} \
                 --pipelineBaseDir ${pipelineBaseDir}