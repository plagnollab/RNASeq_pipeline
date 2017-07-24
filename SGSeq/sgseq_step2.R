## Now take the counts and pretend each event is one gene and each variant is one exon 
## and run in dexseq 

library(SGSeq) 
library(DEXSeq)
library(dplyr) 

library(optparse)
options(echo=T)

option_list <- list(
    make_option(c('--support.tab'), help='', default = ""),
    make_option(c('--step'), help='', default = ""), 
    make_option(c('--code'), help='', default = ""),
    make_option(c('--output.dir'), help='', default = ""),
    make_option(c('--species'), help='', default = "mouse")
)


option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

support.tab <- opt$support.tab
code <- opt$code
step <- opt$step 
output.dir <- opt$output.dir
species <- opt$species 

if(species == "mouse") { 
annotations.tab <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/biomart_annotations_mouse.tab" 
} else if (species == "human") { 
annotations.tab <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/biomart_annotations_human.tab" 
} 


if (step == "step2a") { 
   sgvc.data <- paste0(output.dir, "/", code, "_sgvc.RData") 
   dexseq.data <- paste0(output.dir, "/", code, "_dexseq.RData") 
   res.clean.data <- paste0(output.dir, "/", code, "_res_clean.RData") 
   res.clean.fname <- paste0(output.dir, "/", code, "_res_clean.tab") 
   load(sgvc.data) 
} else if (step == "step2b") { 
   sgvc.data <- paste0(output.dir, "/", code, "_sgv_novel.RData") 
   dexseq.data <- paste0(output.dir, "/", code, "_dexseq_novel.RData") 
   res.clean.data <- paste0(output.dir, "/", code, "_res_clean_novel.RData") 
   res.clean.fname <- paste0(output.dir, "/", code, "_res_clean_novel.tab") 
   load(sgvc.data) 
   sgvc <- sgvc_novel 
} else { 
   message("step needs to be either 2a for known variants or 2b for novel variants") 
}  


varcounts <- counts(sgvc)
vid <- variantID(sgvc)
eid <- eventID(sgvc)

noreads <- which(rowSums(varcounts) == 0 | rowSums(is.na(varcounts)) > 0)

if(length(noreads) > 0) { 
varcounts <- varcounts[-noreads,] 
vid <- vid[-noreads] 
eid <- eid[-noreads] 
} 

print(dim(varcounts)) 
support  <- read.table(support.tab, header = T, stringsAsFactor = F) 

formuladispersion <- ~ sample + (condition + type) * exon
formula0 <-  ~ sample + condition
formula1 <-  ~ sample + condition * exon

if(!file.exists(dexseq.data)) { 
DexSeqExons.loc <- DEXSeqDataSet(countData = varcounts,
                                              sampleData = support,
                                              design = formula1,
                                        featureID = as.factor(vid),
                                        groupID = as.factor(eid))

DexSeqExons.loc <- estimateSizeFactors(DexSeqExons.loc)
DexSeqExons.loc <- DEXSeq::estimateDispersions(DexSeqExons.loc)
DexSeqExons.loc <- DEXSeq::testForDEU(DexSeqExons.loc)
DexSeqExons.loc <- DEXSeq::estimateExonFoldChanges(DexSeqExons.loc)

save(DexSeqExons.loc, file = dexseq.data) 
} else { 
load(dexseq.data) 
} 

res <- DEXSeq::DEXSeqResults (DexSeqExons.loc)
res.clean <- as.data.frame(res)

sample.data <- colData(sgvc)
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

annotation <- read.table(annotations.tab, header = TRUE) 
row.names(annotation) <- annotation$EnsemblID

makePretty <- function(x) { paste(unlist(x), collapse = "+") } 

# Now match the ensembl ID back to more human readable IDs 

getHugoName <- function(x) { paste(annotation[unlist(x), "external_gene_name"], collapse="+" ) } 

res.clean <- dplyr::mutate(res.clean , ensemblName=unlist(lapply(geneName, makePretty))) 
res.clean <- dplyr::mutate(res.clean , geneName=unlist(lapply(geneName, getHugoName))) 
res.clean <- dplyr::mutate(res.clean , txName=unlist(lapply(txName, makePretty))) 
res.clean <- dplyr::mutate(res.clean , variantType=unlist(lapply(variantType, makePretty))) 

save(res.clean, file =  res.clean.data) 
write.table(res.clean, file = res.clean.fname, sep = "\t")
