library(SGSeq) 
library(optparse)
library(stringr)
options(echo=T)
nCores <- 4
option_list <- list(
    make_option(c('--support.tab'), help='', default = "/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i_new_support.tab"),
    make_option(c('--code'), help='', default = "F210I_embryonic_brain_norm"),
    make_option(c('--gtf'), help='', default="/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gtf"),
    make_option(c('--sgseq.anno'), help='', default="/SAN/vyplab/IoN_RNAseq/Kitty/Reference/Mus_musculus.GRCm38.82_sgseq_anno.RData"), 
    make_option(c('--output.dir'), help='', default=""),
    make_option(c('--species'), help='', default="mouse") 
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

support.tab <- opt$support.tab
code <- opt$code
gtf <- opt$gtf
sgseq.anno <- opt$sgseq.anno 
output.dir <- opt$output.dir 
species <- opt$species 

# No longer needed as sgseq.anno is now a required item
# if(species == "mouse") { 
#    sgseq.anno = "/SAN/vyplab/IoN_RNAseq/Kitty/Reference/Mus_musculus.GRCm38.82_sgseq_anno.RData"
# } else if (species == "human") { 
#    sgseq.anno = "/SAN/vyplab/IoN_RNAseq/Kitty/Reference/Homo_sapiens.GRCh38_sgseq_anno.RData"
# } else { 
#    exit("need to put in species") 
# } 
 
# info.file <- paste0(output.dir, "/", code, "_info.RData")
# if(!file.exists(info.file)) { 
# sample.tab <- read.table(support.tab, header = T, stringsAsFactor = F) 
# sample.info <- getBamInfo(sample.tab, cores = nCores) #, yieldSize = 10000 ) 
# save(sample.info, file = paste0(output.dir, "/", code, "_info.RData") ) 
# } else { 
# load(info.file) 
# } 

sample.info.file <- paste0(output.dir, "/", code, "_info.RData")

if(file.exists(sample.info.file)) { 
   load(sample.info.file) 
} else { 
   sample.tab <- read.table(support.tab, header = T, stringsAsFactor = F) 
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




message("loading annotation")
load(sgseq.anno)
print(sgv) 



message("deal with mitochondrial chromosome name")
# chrM in bam files but chrMT in GTF files - stupid!
seqlevels(sgv)[ which(seqlevels(sgv) == "chrMT") ] <- "chrM"


message("getting variant counts") 
print(sample.info) 
sgvc <- getSGVariantCounts(sgv, sample_info = sample.info, cores = nCores, min_denominator = 10, verbose = TRUE)
save(sgvc, file = paste0(output.dir, "/", code, "_sgvc.RData")) 
