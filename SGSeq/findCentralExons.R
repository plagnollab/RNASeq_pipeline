# find central cassette exon coordinates
# from sgf R object produced by SGSeq
library(dplyr)
library(stringr)
library(ggplot2)
library(optparse)
library(SGSeq)
library(data.table)
############
# FUNCTIONS
############

# for each cassette exon
# look at the variant which includes the exon
# this is made up of 3 features
# JEJ - junction exon junction
# the 2nd feature is the exon, which stores the coordinates
findExonPos <- function(varID) { 
  # find the feature IDs for the variant
  featureID <- subset(mcols(sgv_novel), variantID == varID)$featureID     
  # funny regex stuff
  subFeatureID <- gsub("\\(.*\\)", "",featureID) 
  
  if ( subFeatureID == "" ) { 
    subFeatureID = gsub("\\(.*\\),", ",",featureID) 
  }     
  # double check there are only 3 features.
  if( length( str_split(subFeatureID, ",")[[1]] ) != 3){
    message("variant does not contain a central cassete exon")
    return(NULL)
  }
  
  exonID <- strsplit(subFeatureID, ",")[[1]][2]
  
  # Now find the feature that corresponds to this exonID 
  exon_f = sgf_novel[featureID(sgf_novel) == exonID]
  chr = as.character(seqnames(exon_f))  
  start = as.integer(start(exon_f)) 
  end = as.integer(end(exon_f)) 
  strand = as.character(strand(exon_f)) 
  
  pos = c(chr, start, end, strand)  
  
  return(pos) 
} 



# option_list <- list(
#   make_option(c('--support.tab'), help='', default = "/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i_new_support.tab"),
#   make_option(c('--code'), help='', default = "F210I_embryonic_brain_norm"),
#   make_option(c('--case.condition'), help='', default = "HOM"), # now deprecated
#   make_option(c('--sgf_object')),
#   make_option(c('--sgseq_res'), help='', default="/SAN/vyplab/IoN_RNAseq/Kitty/Reference/Mus_musculus.GRCm38.82_sgseq_anno.RData"),
#   make_option(c('--output.dir'), help='', default="")
# )
# 
# option.parser <- OptionParser(option_list=option_list)
# opt <- parse_args(option.parser)
# 
# sgf_object <- opt$sgf_object
# support.tab <- opt$support.tab
# code <- opt$code
# gtf <- opt$gtf
# sgseq.anno <- opt$sgseq.anno 
# output.dir <- opt$output.dir 
# species <- opt$species 
# 

support.tab <- "/Users/Jack/SAN/IoN_RNAseq/Nicol_FUS/mRNAseq/d14/nicol_d14_SGSeq_support.tab"
code <- "Nicol_FUS_d14"
output.dir <- "/Users/Jack/SAN/IoN_RNAseq/Nicol_FUS/mRNAseq/d14/SGSeq/"
step <- "step2b"


support.tab <- "//Users/Jack/SAN/IoN_RNAseq/BilalMalik/SGSeq/12mnth_SGSeq_support.tab"
code <- "Bilal_12mnth"
output.dir <- "//Users/Jack/SAN/IoN_RNAseq/BilalMalik/SGSeq/12mnth/"
step <- "step2b"

###################
## PARSE ARGUMENTS
###################

library(optparse)
options(echo=TRUE)

option_list <- list(
  make_option(c('--support.tab'), help='', default = ""),
  make_option(c('--step'), help='', default = ""), 
  make_option(c('--code'), help='', default = ""),
  make_option(c('--output.dir'), help='', default = ""),
  make_option(c('--species'), help='', default = "")
)


option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

support.tab <- opt$support.tab
code <- opt$code
step <- opt$step 
output.dir <- opt$output.dir
species <- opt$species


support <- read.table(support.tab, header=TRUE, stringsAsFactors = FALSE)
list.conditions <- grep(names(support), pattern = '^condition.*', value  = TRUE)

# find sgv
if (step == "step2a") {
  res.clean.sgf <- paste0(output.dir, "/", code , "_sgv.RData")
  sgv_object <- "sgv"
} else if (step == "step2b") { 
  res.clean.sgf <- paste0(output.dir, "/", code , "_sgv_novel.RData")
  sgv_object <- "sgv_novel"
} else { 
  message("step needs to be either 2a for known variants or 2b for novel variants") 
} 

# find SGF_object
stopifnot( file.exists(res.clean.sgf) )
load(res.clean.sgf)
# doesn't always exist depending on when SGSeq was run
if( !exists(sgv_object) ){
  stop( "sgv cannot be found. You may have to rerun SGSeq as previous versions did not save the sgv object.")
}

for (condition in list.conditions) {
  
  # create a string to represent comparison
  conditions <- support[, condition]
  # remove NA values and find unique
  conditions <- unique(conditions[!is.na(conditions)])
  conditions.name <- paste(conditions, collapse="_")
  
  # sample_list <- support[,c("sample_name",condition)] %>% 
  #   filter(complete.cases(.)) %>% 
  #   split(.,f = .$condition) %>% 
  #   purrr::map("sample_name")
  
  # make condition specific outputs
  condition.dir <- paste0(output.dir,"/", conditions.name )
  
  if (step == "step2a") {
    dexseq.data <- paste0(condition.dir, "/", code, "_", conditions.name, "_dexseq.RData") 
    res.clean.data <- paste0(condition.dir, "/", code, "_", conditions.name,  "_res_clean.RData") 
    res.clean.fname <- paste0(condition.dir, "/", code, "_", conditions.name, "_res_clean.tab")
    res.clean.sgf <- paste0(output.dir, "/", code , "_sgv.RData")
  } else if (step == "step2b") { 
    dexseq.data <- paste0(condition.dir, "/", code, "_", conditions.name, "_dexseq_novel.RData") 
    res.clean.data <- paste0(condition.dir, "/", code, "_", conditions.name,  "_res_clean_novel.RData") 
    res.clean.fname <- paste0(condition.dir, "/", code, "_", conditions.name,  "_res_clean_novel.tab")
    res.clean.sgf <- paste0(output.dir, "/", code , "_sgv_novel.RData")
  } else { 
    message("step needs to be either 2a for known variants or 2b for novel variants") 
  } 
  
  # find SGSeq results
  sgseq_res <- res.clean.fname
  print("reading in SGSeq results")
  #d <- read.table(sgseq_res, header=TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
  d <- fread(sgseq_res, data.table=FALSE)

 
  
  # get all included cassette exons
  cassette_exons <- subset(d, variantType == "SE:I")
  
  # apply function -  very slow
  central_exons <- lapply( cassette_exons$featureID, FUN = findExonPos)
  
  # this removes null entries
  names(central_exons) <- cassette_exons$featureID
  all <- as.data.frame(do.call( rbind, central_exons))
  names(all) <- c("exon.chr", "exon.start", "exon.end", "exon.strand")
  all$exon.start <- as.numeric(as.character(all$exon.start))
  all$exon.end <- as.numeric(as.character(all$exon.end))
  all$featureID <- row.names(all)
  head(all)
  # merge together
  cassette_exons_all <- merge(cassette_exons, all)
  
  # groupID should be kept for further work
  
  # significant cassette exons
  cassette_exons_sig <- filter(cassette_exons_all, padj < 0.05 & !is.na(external_gene_ID) ) %>%
    arrange( padj ) %>%
    select( exon.chr, exon.start, exon.end, external_gene_ID, groupID, exon.strand)
  
  # null cassette exons
  cassette_exons_null <- filter(cassette_exons_all, padj > 0.95 & !is.na(external_gene_ID) ) %>%
    select( exon.chr, exon.start, exon.end, external_gene_ID, groupID, exon.strand)
  
  outFolder <- paste0(condition.dir, "/bed_files")
  
  if( !dir.exists(outFolder)){dir.create(outFolder)}
  
  # write out files
  write.table( cassette_exons_sig, file =  paste0(outFolder, "/cassette_exons_central_sig.bed"), col.names=FALSE, row.names=FALSE, quote = FALSE , sep = "\t")
  write.table( cassette_exons_null, file =  paste0(outFolder, "/cassette_exons_central_null.bed"), col.names=FALSE, row.names=FALSE, quote = FALSE , sep = "\t")
  
}
