# Builds a list of variants for sgseq from a gff file 
library(SGSeq) 
library(optparse)

option_list <- list(
    make_option(c('--gtf'), help='', default="/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gtf"),
    make_option(c('--sgseq.anno'), help='', default="/SAN/vyplab/IoN_RNAseq/Kitty/Reference/Mus_musculus.GRCm38.82_sgseq_anno.RData")#, 
    #make_option(c('--refFolder'), help='', default=""),
#    make_option(c('--species'), help='', default="mouse") 
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

gtf <- opt$gtf
sgseq.anno <- opt$sgseq.anno 
#output.dir <- opt$output.dir 
#species <- opt$species 

# create outFile from output dir and GTF name
# deprecated - use sgseq.anno flag instead
#outFile <- paste0( refFolder, gsub("gtf", "SGSeq_anno.Rdata", basename(gtf) ) )

print( "creating SGSeq transcript annotations")
print( paste( "from", gtf ) )


# tx <- importTranscripts(gtf) 
# txf_ucsc <- convertToTxFeatures(tx)
# sgf_ucsc <- convertToSGFeatures(txf_ucsc)
# sgv <- findSGVariants(sgf_ucsc)

tx <- importTranscripts(gtf) 
txf <- convertToTxFeatures(tx) 
sgf <- convertToSGFeatures(txf) 
sgv <- findSGVariants(sgf)

print( paste( "saving to:", sgseq.anno )  )
save(tx, txf, sgf, sgv, file = sgseq.anno) 
