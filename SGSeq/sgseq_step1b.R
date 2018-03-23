library(SGSeq) 
library(optparse)
library(stringr)
options(echo=T)

# for debugging
#load("../../Fly_C9/SGSeq/Fly_C9_txf_novel.RData")

nCores <- 4
# Finds novel transcripts and obtain the variant counts 

option_list <- list(
    make_option(c('--support.tab'), help='', default = "/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i_new_support.tab"),
    make_option(c('--code'), help='', default = "F210I_embryonic_brain_norm"),
    #make_option(c('--case.condition'), help='', default = "HOM"), # now deprecated
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
   sample.tab <- read.table(support.tab, header = TRUE, stringsAsFactor = FALSE) 
   # this takes forever, give more cores and only look at a subset 
   sample.info <- getBamInfo(sample.tab, cores = nCores , yieldSize = 10000 ) 
   

   logFiles <- gsub("_unique.bam", "Log.final.out", sample.tab$file_bam )

  getLibSize <- function(logFile){
    stopifnot(file.exists(logFile))
    log <- readLines(logFile)
    unique <- log[ grepl("Uniquely mapped reads number", log)]
    multi <- log[ grepl("Number of reads mapped to multiple loci", log)]

    num_unique <- str_trim( str_split_fixed( unique, "\\|\t", 2)[,2] )
    num_multi <- str_trim( str_split_fixed( multi, "\\|\t", 2)[,2] )
    libSize <- as.numeric(num_unique) + as.numeric(num_multi)
    return(libSize)
  }
  sample.info$lib_size <- sapply(logFiles, FUN = getLibSize)


   save(sample.info, file = sample.info.file) 
} 

#tx <- importTranscripts(gtf) 
#txf <- convertToTxFeatures(tx) 
#sgf <- convertToSGFeatures(txf) 
message("loading annotation")
load(sgseq.anno)

#si_cases <- subset(sample.info, condition == case.condition)
#print(si_cases) 
#txf_novel <- predictTxFeatures(si_cases, min_junction_count = 5, verbose = TRUE, cores = nCores)

#####
# this step predicts features for each sample and merges features together
# it's very slow - may be better to do this individually for each sample and merge back together
#####

# originally just did this for the cases and not the controls - not sure about this.
txf_novel <- predictTxFeatures(sample.info, min_junction_count = 5, verbose = FALSE, cores = nCores, max_complexity = 15)

save(txf_novel, file = paste0(output.dir, "/", code, "_txf_novel.RData")) 

# drop any sequence features that aren't part of standard chromosomes
# only works for humans - what about flies? remove anything that lacks "_".
#txf_novel <- keepSeqlevels(txf_novel, paste0("chr",c(1:19, "X", "Y")))
keepChrs <- seqnames(txf_novel@seqinfo)[!grepl( "_", seqnames(txf_novel@seqinfo) ) ]
print( "keeping only sequences from:")
print(keepChrs)
txf_novel <- keepSeqlevels( txf_novel, keepChrs )


sgf_novel <- convertToSGFeatures(txf_novel)
sgf_novel <- annotate(sgf_novel, txf)
sgv_novel <- findSGVariants(sgf_novel, include = "all", cores = nCores)

sgvc_novel <- getSGVariantCounts(sgv_novel, sample_info = sample.info, cores = nCores, min_denominator = 10, verbose = TRUE)

save(sgf_novel, sgv_novel, sgvc_novel, file = paste0(output.dir, "/", code, "_sgv_novel.RData")) 
