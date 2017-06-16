## Now take the counts from polyA site counting 
## and run in dexseq 

library(DEXSeq)
library(dplyr) 
library(BiocParallel)
library(optparse)

options(echo=TRUE)

option_list <- list(
    make_option(c('--support.tab'), help='', default = ""),
    make_option(c('--code'), help='', default = ""),
    make_option(c('--output.dir'), help='', default = ""),
    make_option(c('--input.dir'), help='', default = ""),
    make_option(c('--biomartAnnotation'), help='', default = "mouse"),
    make_option(c('--nCores'), help='', default = 4)
)


option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

support.tab <- opt$support.tab
code <- opt$code
output.dir <- opt$output.dir
input.dir <- opt$input.dir 
biomartAnnotation <- opt$biomartAnnotation
nCores <- opt$nCores

#support.tab="/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/threeprimeseq/d14/d14_hom_wt_support.tab"
#code <- "threeprime_d14_hom_wt" 
#output.dir=paste0("/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/threeprimeseq/d14/dexseq/", code) 
#input.dir="/SAN/vyplab/IoN_RNAseq/Kitty/Nicol/threeprimeseq/d14/"
#species <- "mouse" 

dir.create(output.dir) 
BPPARAM = MulticoreParam(workers=8)

support <- read.table(support.tab, sep ='\t', stringsAsFactors=FALSE, header = TRUE, comment.char = "") 
samples <- support$sample 

# if(species == "mouse") { 
# annotations.tab <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/biomart_annotations_mouse.tab" 
# } else if (species == "human") { 
# annotations.tab <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/biomart_annotations_human.tab" 
# } 

annotations <- read.table(biomartAnnotation, sep = "\t", header = TRUE) 

# Read in each counts file and make a matrix 
for (sample in samples) { 
   bedfile = paste0(output.dir, "/", sample, "/", sample, "_clusters_counts.bed") 
   beddata <- read.table(bedfile, sep = "\t", header = FALSE) 
   counts <- beddata[,8]
   genes <- beddata[,7] 
   message(sample) 
   if (!exists("all.counts")) { 
       all.counts <- as.data.frame(counts)  
       names(all.counts)[length(all.counts)] <- sample
	# need to fix - dexseq doesn't like colons in rownames 
       row.names(all.counts) <- paste0(beddata[,7], ":", beddata[,1], ":", beddata[,2], "-", beddata[,3]) 
   } else {
       all.counts <- cbind(all.counts, as.data.frame(counts) ) 
       names(all.counts)[length(all.counts)] <- sample 
   } 
   print(head(all.counts))    
} 

formuladispersion <- ~ sample + (condition + type) * exon
formula0 <-  ~ sample + condition
formula1 <-  ~ sample + condition * exon

DexSeqExons.loc <- DEXSeqDataSet(countData = all.counts,
                                              sampleData = support,
                                              design = formula1,
                                        featureID = as.factor(row.names(all.counts)),
                                        groupID = as.factor(genes))

DexSeqExons.loc <- estimateSizeFactors(DexSeqExons.loc)
DexSeqExons.loc <- DEXSeq::estimateDispersions(DexSeqExons.loc, BPPARAM=BPPARAM)
DexSeqExons.loc <- DEXSeq::testForDEU(DexSeqExons.loc, BPPARAM=BPPARAM)
DexSeqExons.loc <- DEXSeq::estimateExonFoldChanges(DexSeqExons.loc, BPPARAM=BPPARAM)

dexseq.data <- paste0(output.dir, "/", "dexseq.RData") 
save(DexSeqExons.loc, file = dexseq.data) 

res <- DEXSeq::DEXSeqResults (DexSeqExons.loc)
res <- as.data.frame(res) 

logname <- grep(names(res), pattern = 'log2fold', value = TRUE)
res.clean <- as(res[, c('groupID', 'exonBaseMean', logname, 'dispersion', 'stat', 'pvalue')], 'data.frame')
names(res.clean)<- c("EnsemblID", "meanBase", logname, "dispersion", "stat", "pvalue")

res.clean$FDR <- p.adjust(res.clean$pvalue, method = 'fdr')
res.clean$ID <- row.names(res.clean) 

res.clean$external_gene_name <- annotations$external_gene_name[ match(res.clean$EnsemblID, table = annotations$EnsemblID) ]
res.clean <- res.clean[, c('external_gene_name', "EnsemblID", "meanBase", logname, "dispersion", "stat", "pvalue", "FDR", "ID")]

res.clean <- res.clean[order(res.clean$FDR), ] 

write.table(res.clean, row.names = FALSE, quote = FALSE, sep = "\t", file = paste0(output.dir, "/dexseq_res_clean.tab"))  

