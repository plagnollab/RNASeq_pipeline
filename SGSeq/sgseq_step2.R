## Now take the counts and pretend each event is one gene and each variant is one exon 
## and run in dexseq 
# Kitty Lo wrote this
# Jack Humphrey tweaked it

library(SGSeq) 
library(DEXSeq)
library(dplyr)
library(ggplot2)
library(stringr)

BPPARAM <- MulticoreParam(workers=8)


library(optparse)
options(echo=TRUE)

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

list.covariates <- grep(names(support), pattern = '^type.*', value = TRUE )

if( length(list.covariates) == 1 ){
	add.covariate <- TRUE

}else{
	add.covariate <- FALSE
}

# pick which counts file to use
if (step == "step2a") { 
   sgvc.data <- paste0(output.dir, "/", code, "_sgvc.RData") 
   print(sgvc.data)
   load(sgvc.data) 
} else if (step == "step2b") { 
   sgvc.data <- paste0(output.dir, "/", code, "_sgv_novel.RData") 
   print(sgvc.data)
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

  
  support.loc$condition <- factor(support[, condition])
	  
  # subset the samples from the matrix that are to be tested
  varcounts.loc <- varcounts[, !is.na(support[,condition]) ]
  message("varcounts subset OK")
  #loc.countFiles <- countFiles[ !is.na(support.loc$condition) ]
  support.loc <-  support.loc[ !is.na(support.loc$condition), ]
  message("support subset OK")
  # create formula
  #formuladispersion <- ~ sample + (condition + type) * exon
  #formula0 <-  ~ sample + condition
  #formula1 <-  ~ sample + condition * exon
 
  
 
  if( add.covariate == TRUE){
    names(support.loc)[ which( names(support.loc) == list.covariates) ] <- "type" 
    formulaFullModel = ~ sample + exon + type:exon + condition:exon
    formulaReducedModel = ~ sample + exon + type:exon
  }else{
    formulaFullModel = ~ sample + exon + condition:exon
    formulaReducedModel = ~ sample + exon
  }

# save local environment for testing
message("saving")
    save.image(file = "test.Rdata")	
 

  # run DEXSeq
  if(!file.exists(dexseq.data)) { 
    DexSeqExons.loc <- DEXSeqDataSet(countData = varcounts.loc,
                                     sampleData = support.loc,
                                     #design = formula1,
                                     design = formulaFullModel,
				     featureID = as.factor(vid),
                                     groupID = as.factor(eid))
   
 
    DexSeqExons.loc <- estimateSizeFactors(DexSeqExons.loc)
    DexSeqExons.loc <- DEXSeq::estimateDispersions(DexSeqExons.loc, formula = formulaFullModel, BPPARAM=BPPARAM)
    DexSeqExons.loc <- DEXSeq::testForDEU(DexSeqExons.loc,
			reducedModel = formulaReducedModel,
                  	fullModel = formulaFullModel,
 			BPPARAM=BPPARAM)

    DexSeqExons.loc <- DEXSeq::estimateExonFoldChanges(DexSeqExons.loc, BPPARAM=BPPARAM, fitExpToVar="condition")
    
    # save data
    save(DexSeqExons.loc, sgvc, support, annotation, file = dexseq.data) 
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
  
  # find all columns that correspond to "countData.x"
  # rename them to the naming scheme
  countData.cols <- which( grepl("countData", names(res.clean) ) )
  names(res.clean)[countData.cols] <- sample.names
  #count.start <- which(names(res.clean) == "countData.1") 
  #last <- count.start+n.samples-1 
  #names(res.clean)[count.start:last] <- sample.names 
  
  res.clean$genomicData <- NULL 
  res.clean$FDR <- p.adjust(res.clean$pvalue, method = 'fdr')
  
  sgvc.df <- as.data.frame(mcols(sgvc) )
  psi <- variantFreq(sgvc)
  
  if (length(noreads) > 0) { 
    sgvc.df <- sgvc.df[-noreads,]
    psi <- psi[-noreads,] 
  } 

  # select just the samples present in the condition
  psi <- psi[, sample.names]
  
  colnames(psi) <- paste0(colnames(psi)  , "_psi") 
  res.clean <- cbind(res.clean, psi) 
  res.clean$geneName <- sgvc.df$geneName
  res.clean$txName <- sgvc.df$txName
  res.clean$eventID <- sgvc.df$eventID
  res.clean$variantType <- sgvc.df$variantType
  res.clean$from <- sgvc.df$from
  res.clean$to <- sgvc.df$to
  res.clean$type <- sgvc.df$type

  res.clean$external_gene_ID <- annotation$external_gene_name[ match( res.clean$geneName, annotation$EnsemblID)]

  res.clean <- res.clean[order(res.clean$pvalue), ]

  createCoords <- function(from, to){
    fromSplit <- str_split_fixed(from, ":", 4)
    toSplit <- str_split_fixed(to, ":", 4)
    coord <- ifelse( toSplit[,4] == "+", 
      yes = paste0( toSplit[,2], ":", fromSplit[,3], "-", toSplit[,3]),
      no = paste0( toSplit[,2], ":", toSplit[,3], "-", fromSplit[,3]) 
      ) 
  }
  
  res.clean$coords <- createCoords( from = res.clean$from, to = res.clean$to)

  makePretty <- function(x) { paste(unlist(x), collapse = "+") } 
  
  # Now match the ensembl ID back to more human readable IDs 
  
  getHugoName <- function(x) { paste(annotation[unlist(x), "external_gene_name"], collapse="+" ) } 
  
  res.clean <- dplyr::mutate(res.clean , ensemblName=unlist(lapply(geneName, makePretty))) 
  res.clean <- dplyr::mutate(res.clean , geneName=unlist(lapply(geneName, getHugoName))) 
  res.clean <- dplyr::mutate(res.clean , txName=unlist(lapply(txName, makePretty))) 
  res.clean <- dplyr::mutate(res.clean , variantType=unlist(lapply(variantType, makePretty))) 
  
  # save and write
  save(res.clean, file =  res.clean.data) 
  write.table(res.clean, file = res.clean.fname, sep = "\t", row.names = FALSE)
  
  source("/SAN/vyplab/HuRNASeq/RNASeq_pipeline/SGSeq/makePieChartsAllEvents.R")
  
  try(
    makePieChart(sgseqRes = res.clean, title = paste0(code, "_", conditions.name), FDRlimit = 0.05, outFolder = condition.dir )
  )
}

