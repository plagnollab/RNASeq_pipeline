########## first the high level parameters about where the scripts and the bundle are located
RNASEQPIPBASE=/SAN/vyplab/HuRNASeq/RNASeq_pipeline/ 
RNASEQBUNDLE=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle  ##this probably should stay

export RNASEQPIPBASE=$RNASEQPIPBASE
export RNASEQBUNDLE=$RNASEQBUNDLE
pipeline=${RNASEQPIPBASE}/RNAseq_pipeline_v8.sh

#############
species=human_hg38
submit=yes
force=no
step0_QC=no

trim_galore=no
starStep1a=no
starStep1b=no
starStep2=yes

summary=no
prepareCounts=yes
Rdeseq=yes
Rdexseq=yes
goseq=yes

stranded=fr-secondstrand

oFolder=/SAN/vyplab/HuRNASeq/ENCODE/SFPQ/K562_ENCSR535YPK/processed
iFolder=/SAN/vyplab/HuRNASeq/ENCODE/
dataframe=/SAN/vyplab/HuRNASeq/ENCODE/SFPQ/K562_ENCSR535YPK/K562_ENCSR535YPK_support.tab
code=SFPQ_K562_ENCSR535YPK

sh -x $pipeline --goseq ${goseq}  --step0_QC $step0_QC --trim_galore $trim_galore --iFolder ${iFolder} --oFolder ${oFolder} --dataframe ${dataframe} --code ${code} --prepareCounts ${prepareCounts} --Rdexseq ${Rdexseq} --Rdeseq ${Rdeseq}  --summary ${summary} --starStep1a ${starStep1a} --starStep1b ${starStep1b} --starStep2 ${starStep2} --species ${species} --submit ${submit} --force ${force}




echo $mainscript







