## Now take the counts and pretend each event is one gene and each variant is one exon 
## and run in dexseq 
# Kitty Lo wrote this
# Jack Humphrey tweaked it

library(SGSeq) 
library(DEXSeq)
library(dplyr)
library(ggplot2)

library(optparse)
options(echo=T)

option_list <- list(
    make_option(c('--support.tab'), help='', default = ""),
    make_option(c('--step'), help='', default = ""), 
    make_option(c('--code'), help='', default = ""),
    make_option(c('--output.dir'), help='', default = ""),
    make_option(c('--annotation'), help='', default = "")
)


option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

support.tab <- opt$support.tab
code <- opt$code
step <- opt$step 
output.dir <- opt$output.dir
species <- opt$species
annotations.tab <- opt$annotation

# if(species == "mouse") { 
# annotations.tab <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/biomart_annotations_mouse.tab" 
# } else if (species == "human") { 
# annotations.tab <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/biomart_annotations_human.tab" 
# } 

# read in annotation
annotation <- read.table(annotations.tab, header = TRUE) 
row.names(annotation) <- annotation$EnsemblID


# read in support
support  <- read.table(support.tab, header = T, stringsAsFactor = F) 
list.conditions <- grep(names(support), pattern = '^condition.*', value  = TRUE)

# pick which counts file to use
if (step == "step2a") { 
   sgvc.data <- paste0(output.dir, "/", code, "_sgvc.RData") 
   load(sgvc.data) 
} else if (step == "step2b") { 
   sgvc.data <- paste0(output.dir, "/", code, "_sgv_novel.RData") 
   load(sgvc.data) 
   sgvc <- sgvc_novel 
} else { 
   message("step needs to be either 2a for known variants or 2b for novel variants") 
}  

# this is the counts matrix
varcounts <- counts(sgvc)
vid <- variantID(sgvc)
eid <- eventID(sgvc)

# for testing
##varcounts <- varcounts[1:100,]
##vid <- vid[1:100]
##eid <- eid[1:100]



# remove no count rows
noreads <- which(rowSums(varcounts) == 0 | rowSums(is.na(varcounts)) > 0)
if(length(noreads) > 0) { 
varcounts <- varcounts[-noreads,] 
vid <- vid[-noreads] 
eid <- eid[-noreads] 
} 

print(dim(varcounts)) 

for (condition in list.conditions) {
  
  # create a string to represent comparison
  conditions <- support[, condition]
  # remove NA values and find unique
  conditions <- unique(conditions[!is.na(conditions)])
  conditions.name <- paste(conditions, collapse="_")
   
  # make condition specific outputs
  condition.dir <- paste0(output.dir,"/", conditions.name )
  
  if( !dir.exists(condition.dir) ){ 
      dir.create(condition.dir) 
  }
     
  if (step == "step2a") {
   dexseq.data <- paste0(condition.dir, "/", code, "_", conditions.name, "_dexseq.RData") 
   res.clean.data <- paste0(condition.dir, "/", code, "_", conditions.name,  "_res_clean.RData") 
   res.clean.fname <- paste0(condition.dir, "/", code, "_", conditions.name, "_res_clean.tab") 
  } else if (step == "step2b") { 
       dexseq.data <- paste0(condition.dir, "/", code, "_", conditions.name, "_dexseq_novel.RData") 
       res.clean.data <- paste0(condition.dir, "/", code, "_", conditions.name,  "_res_clean_novel.RData") 
       res.clean.fname <- paste0(condition.dir, "/", code, "_", conditions.name,  "_res_clean_novel.tab")  
  } else { 
   message("step needs to be either 2a for known variants or 2b for novel variants") 
  }  

  message('Condition ', condition)
  support.loc <- support

  ##handle the type variable
  support.loc$condition <- factor(support[, condition])
  
  # subset the samples from the matrix that are to be tested
  varcounts.loc <- varcounts[, !is.na(support[,condition]) ]
  #loc.countFiles <- countFiles[ !is.na(support.loc$condition) ]
  support.loc <-  support.loc[ !is.na(support.loc$condition), ]

  # create formula
  formuladispersion <- ~ sample + (condition + type) * exon
  formula0 <-  ~ sample + condition
  formula1 <-  ~ sample + condition * exon
  
  # run DEXSeq
  if(!file.exists(dexseq.data)) { 
    DexSeqExons.loc <- DEXSeqDataSet(countData = varcounts.loc,
                                     sampleData = support.loc,
                                     design = formula1,
                                     featureID = as.factor(vid),
                                     groupID = as.factor(eid))
    
    DexSeqExons.loc <- estimateSizeFactors(DexSeqExons.loc)
    DexSeqExons.loc <- DEXSeq::estimateDispersions(DexSeqExons.loc)
    DexSeqExons.loc <- DEXSeq::testForDEU(DexSeqExons.loc)
    DexSeqExons.loc <- DEXSeq::estimateExonFoldChanges(DexSeqExons.loc)
    
    # save data
    save(DexSeqExons.loc, file = dexseq.data) 
  } else { # mainly for testing
    load(dexseq.data) 
  } 
  
  res <- DEXSeq::DEXSeqResults (DexSeqExons.loc)
  res.clean <- as.data.frame(res)
  
  ## just for testing
  ##}
  ##quit()
  
  sample.data <- colData(sgvc)
  # remove any not in our condition
  sample.data <- sample.data[!is.na(support[,condition]), ]
  
  sample.names <- sample.data$sample_name 
  n.samples <- length(sample.names) 
  count.start <- which(names(res.clean) == "countData.1") 
  last <- count.start+n.samples-1 
  names(res.clean)[count.start:last] <- sample.names 
  
  res.clean$genomicData <- NULL 
  res.clean$FDR <- p.adjust(res.clean$pvalue, method = 'fdr')
  
  sgvc.df <- as.data.frame(mcols(sgvc) )
  psi <- variantFreq(sgvc)
  
  if (length(noreads) > 0) { 
    sgvc.df <- sgvc.df[-noreads,]
    psi <- psi[-noreads,] 
  } 
  
  colnames(psi) <- paste0(colnames(psi)  , "_psi") 
  res.clean <- cbind(res.clean, psi) 
  res.clean$geneName <- sgvc.df$geneName
  res.clean$txName <- sgvc.df$txName
  res.clean$eventID <- sgvc.df$eventID
  res.clean$variantType <- sgvc.df$variantType
  res.clean$from <- sgvc.df$from
  res.clean$to <- sgvc.df$to
  res.clean$type <- sgvc.df$type
  res.clean <- res.clean[order(res.clean$pvalue), ]
  
  makePretty <- function(x) { paste(unlist(x), collapse = "+") } 
  
  # Now match the ensembl ID back to more human readable IDs 
  
  getHugoName <- function(x) { paste(annotation[unlist(x), "external_gene_name"], collapse="+" ) } 
  
  res.clean <- dplyr::mutate(res.clean , ensemblName=unlist(lapply(geneName, makePretty))) 
  res.clean <- dplyr::mutate(res.clean , geneName=unlist(lapply(geneName, getHugoName))) 
  res.clean <- dplyr::mutate(res.clean , txName=unlist(lapply(txName, makePretty))) 
  res.clean <- dplyr::mutate(res.clean , variantType=unlist(lapply(variantType, makePretty))) 
  
  # save and write
  save(res.clean, file =  res.clean.data) 
  write.table(res.clean, file = res.clean.fname, sep = "\t")
  
  source("/SAN/vyplab/HuRNASeq/RNASeq_pipeline/SGSeq/makePieChartsAllEvents.R")
  
  try(
    makePieChart(sgseqRes = res.clean, title = paste0(code, "_", conditions.name), FDRlimit = 0.01, outFolder = condition.dir )
  )
}

