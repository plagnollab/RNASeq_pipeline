pipelineBaseDir=/cluster/project8/vyp/kitty/RNASeq_pipeline/juncseq/
qorts_script=${pipelineBaseDir}/run_qorts.sh

gtf=/SAN/vyplab/IoN_RNAseq/Kitty/Reference/gencode.vM11.chr_patch_hapl_scaff.basic.annotation.gtf
basedir=/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/ko/
support=/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/ko/nicol_fus_ko_support.tab
baseoutputdir=/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/ko/
code=nicol_ko 
submit=yes

if [ ! -e $outputDir ]; then
                mkdir $outputDir
fi

sh $qorts_script --gtf ${gtf} \
                 --support ${support} \
                 --submit ${submit} \
                 --baseoutputdir ${baseoutputdir} \
                 --basedir ${basedir} 
