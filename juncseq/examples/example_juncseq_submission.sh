pipelineBaseDir=/cluster/project8/vyp/kitty/RNASeq_pipeline/juncseq/
juncseq_script=${pipelineBaseDir}/run_juncseq.sh
juncseq_rscript=${pipelineBaseDir}/juncseq.R

gtf=/SAN/vyplab/IoN_RNAseq/Kitty/Reference/gencode.vM11.chr_patch_hapl_scaff.basic.annotation.gtf
gff=/SAN/vyplab/IoN_RNAseq/Kitty/Reference/gencode.vM11.chr_patch_hapl_scaff.basic.annotation.juncseq.gff
basedir=/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/ko/
support=/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/ko/nicol_fus_ko_support.tab
baseoutputdir=/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/ko/
code=nicol_ko 
submit=yes

if [ ! -e $outputDir ]; then
                mkdir $outputDir
fi

sh $juncseq_script --gff ${gff} \
                 --support ${support} \
                 --submit ${submit} \
                 --baseoutputdir ${baseoutputdir} \
                 --basedir ${basedir} \
                 --code ${code} \
                 --juncseqscript ${juncseq_rscript} 
