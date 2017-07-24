# Builds a list of variants for sgseq from a gff file 
library(SGSeq) 

tx <- importTranscripts("/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/Homo_sapiens.GRCh38.82_fixed.gtf") 
txf_ucsc <- convertToTxFeatures(tx)
sgf_ucsc <- convertToSGFeatures(txf_ucsc)
sgv <- findSGVariants(sgf_ucsc)
save(sgv, file = "/SAN/vyplab/IoN_RNAseq/Kitty/Reference/Homo_sapiens.GRCh38_sgseq_anno.RData") 
