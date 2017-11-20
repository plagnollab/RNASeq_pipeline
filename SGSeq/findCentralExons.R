# find central cassette exon coordinates
# from sgf R object produced by SGSeq
library(dplyr)
library(stringr)
library(ggplot2)
library(optparse)

option_list <- list(
  make_option(c('--support.tab'), help='', default = "/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i_new_support.tab"),
  make_option(c('--code'), help='', default = "F210I_embryonic_brain_norm"),
  make_option(c('--case.condition'), help='', default = "HOM"), # now deprecated
  make_option(c('--sgf_object')),
  make_option(c('--sgseq_res'), help='', default="/SAN/vyplab/IoN_RNAseq/Kitty/Reference/Mus_musculus.GRCm38.82_sgseq_anno.RData"),
  make_option(c('--output.dir'), help='', default="")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

sgf_object <- opt$sgf_object
support.tab <- opt$support.tab
code <- opt$code
gtf <- opt$gtf
sgseq.anno <- opt$sgseq.anno 
output.dir <- opt$output.dir 
species <- opt$species 



load(sgf_object)
d <- as.data.frame(fread(sgseq_res), stringsAsFactors=FALSE)

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

# get all included cassette exons
cassette_exons <- subset(d, variantType == "SE:I")
# very slow
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

outFolder <- paste0(dirname(sgseq_res), "/bed_files")
if( !dir.exists(outFolder)){dir.create(outFolder)}

# write out files
write.table( cassette_exons_sig, file =  paste0(outFolder, "/cassette_exons_central_sig.bed"), col.names=FALSE, row.names=FALSE, quote = FALSE , sep = "\t")
write.table( cassette_exons_null, file =  paste0(outFolder, "/cassette_exons_central_null.bed"), col.names=FALSE, row.names=FALSE, quote = FALSE , sep = "\t")

