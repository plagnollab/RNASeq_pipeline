########## first the high level parameters about where the scripts and the bundle are located
RNASEQPIPBASE=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline  ##this probably needs to be updated, but I do not think it has to
RNASEQBUNDLE=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle  ##this probably should stay

export RNASEQPIPBASE=$RNASEQPIPBASE
export RNASEQBUNDLE=$RNASEQBUNDLE
pipeline=${RNASEQPIPBASE}/RNAseq_pipeline_v8.sh

#############
species=Human_hg38
submit=yes
force=yes

trim_galore=yes

step0_QC=yes
starStep1a=yes
starStep1b=yes
starStep2=yes

summary=no
prepareCounts=no
Rdeseq=no
Rdexseq=no



oFolder=/scratch2/vyp-scratch2/IoO_RNASeq/processed/Owen_blood
iFolder=/SAN/biosciences/ioo_no_aw/data/Project/A2914/fastq

dataframe=support/Owen_blood.tab
code=Owen_blood
species=human_hg38

sh $pipeline --iFolder ${iFolder} --oFolder ${oFolder} --dataframe ${dataframe} --code ${code} --prepareCounts ${prepareCounts} --Rdexseq ${Rdexseq} --Rdeseq ${Rdeseq}  --summary ${summary} --step0_QC ${step0_QC} --trim_galore ${trim_galore} --starStep1a ${starStep1a} --starStep1b ${starStep1b} --starStep2 ${starStep2} --species ${species} --submit ${submit} --force ${force}




echo $mainscript







