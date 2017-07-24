library(SGSeq) 
library(optparse)
options(echo=T)

# Finds novel transcripts and obtain the variant counts 

option_list <- list(
    make_option(c('--support.tab'), help='', default = "/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i_new_support.tab"),
    make_option(c('--code'), help='', default = "F210I_embryonic_brain_norm"),
    make_option(c('--case.condition'), help='', default = "HOM"),
    make_option(c('--sgseq.anno'), help='', default="/SAN/vyplab/IoN_RNAseq/Kitty/Reference/Mus_musculus.GRCm38.82_sgseq_anno.RData"),
    make_option(c('--gtf'), help='', default="/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gtf"),
    make_option(c('--output.dir'), help='', default="")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

support.tab <- opt$support.tab
code <- opt$code
gtf <- opt$gtf
case.condition <- opt$case.condition
sgseq.anno <- opt$sgseq.anno 
output.dir <- opt$output.dir 
sample.info.file <- paste0(output.dir, "/", code, "_info.RData")

if(file.exists(sample.info.file)) { 
   load(sample.info.file) 
} else { 
   sample.tab <- read.table(support.tab, header = T, stringsAsFactor = F) 
   # this takes forever, give more cores and only look at a subset 
   sample.info <- getBamInfo(sample.tab, cores = 4, yieldSize = 10000 ) 
   save(sample.info, file = sample.info.file) 
} 

#tx <- importTranscripts(gtf) 
#txf <- convertToTxFeatures(tx) 
#sgf <- convertToSGFeatures(txf) 
message("loading annotation")
load(sgseq.anno)

si_cases <- subset(sample.info, condition == case.condition)
print(si_cases) 
txf_novel <- predictTxFeatures(si_cases, min_junction_count = 5, verbose = TRUE, cores = 4)

save(txf_novel, file = paste0(output.dir, "/", code, "_txf_novel.RData")) 

txf_novel <- keepSeqlevels(txf_novel, paste0("chr",c(1:19, "X", "Y")))
sgf_novel <- convertToSGFeatures(txf_novel)
sgf_novel <- annotate(sgf_novel, txf)
sgv_novel <- findSGVariants(sgf_novel, include = "all")

sgvc_novel <- getSGVariantCounts(sgv_novel, sample_info = sample.info, cores = 8, min_denominator = 10, verbose = TRUE)
save(sgf_novel, sgv_novel, sgvc_novel, file = paste0(output.dir, "/", code, "_sgv_novel.RData")) 
