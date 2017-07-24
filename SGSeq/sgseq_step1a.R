library(SGSeq) 
library(optparse)
options(echo=T)

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

if(species == "mouse") { 
   sgseq.anno = "/SAN/vyplab/IoN_RNAseq/Kitty/Reference/Mus_musculus.GRCm38.82_sgseq_anno.RData"
} else if (species == "human") { 
   sgseq.anno = "/SAN/vyplab/IoN_RNAseq/Kitty/Reference/Homo_sapiens.GRCh38_sgseq_anno.RData"
} else { 
   exit("need to put in species") 
} 
 
info.file <- paste0(output.dir, "/", code, "_info.RData")
if(!file.exists(info.file)) { 
sample.tab <- read.table(support.tab, header = T, stringsAsFactor = F) 
sample.info <- getBamInfo(sample.tab) 
save(sample.info, file = paste0(output.dir, "/", code, "_info.RData") ) 
} else { 
load(info.file) 
} 

message("loading annotation")
load(sgseq.anno)
print(sgv) 

message("getting variant counts") 
print(sample.info) 
sgvc <- getSGVariantCounts(sgv, sample_info = sample.info, cores = 4, min_denominator = 10, verbose = TRUE)
save(sgvc, file = paste0(output.dir, "/", code, "_sgvc.RData")) 
