library(dplyr)
library(stringr)
library(data.table)
library(optparse)

options(echo=TRUE)

option_list <- list(
    make_option(c('--code'), help='', default = ""),
    make_option(c('--output.dir'), help='', default = ""),
    make_option(c('--biomartAnnotation'), help='', default = "mouse"),
    make_option(c('--dexseqResults'), help='', default = "mouse"),
    make_option(c('--mode'), help='')
)


option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

code <- opt$code
output.dir <- opt$output.dir
biomartAnnotations <- opt$biomartAnnotation
dexseq.results <- opt$dexseqResults
mode <- opt$mode



# iFolder <- "/Users/Jack/google_drive/TDP_paper/"
# biomartAnnotations <- paste0(iFolder, "/differential_splicing/biomart_annotations_mouse.tab")
biomart <- as.data.frame(fread(biomartAnnotations))
#dexseq.results <- "/Users/Jack/SAN/IoN_RNAseq/Nicol/polyA/FUS_ko_HOM/results/dexseq_results_FUS_ko_HOM_thinned.tab"

res <- as.data.frame(fread(dexseq.results))

# get the clusters for all genes that have at least one significant cluster
sigClusters <- filter(res, FDR < 0.1)
sigGenes <- unique(sigClusters$external_gene_name, na.rm=TRUE)
sigGenes <- sigGenes[ !is.na(sigGenes)]


allClusters <- filter(res, external_gene_name %in% sigGenes )
coords <- str_split_fixed( str_split_fixed(allClusters$ID, ":", 3)[,3], ":", 2)
coords <- data.frame(
  chr = coords[,1],
  start = str_split_fixed(coords[,2], "-", 2)
)
names(coords) <- c("chr","start","end")

allClusters <- cbind(allClusters,coords)
allClusters$strand <- biomart$strand[ match(allClusters$EnsemblID, biomart$EnsemblID)]

names(allClusters)[4] <- "log2FoldChange"

# for each gene
#for( i in 1:length(sigGenes) ){

clusterDirections <- lapply(1:length(sigGenes), FUN = function(i){ 
  cluster <- filter(allClusters, external_gene_name == sigGenes[i])
  # assign proximal and distal labels
  # if only 2 this is easy
  n <- nrow(cluster)
  #print(n)
  # order by position so proximal clusters at top
  if( unique(cluster$strand) == 1){
  cluster <- cluster[order(cluster$start),]
  }else{
    cluster <- cluster[order(cluster$start, decreasing = TRUE),]
  }
  cluster$label <- NA 
  cluster$label[1] <- "proximal"
  cluster$label[n] <- "distal"
  if( n > 3){
    midpoint <- floor(n/2)
    cluster$label[1:midpoint] <- "proximal"
    cluster$label[(midpoint+1):n ] <- "distal"
  }
  sigCluster <- subset(cluster, FDR < 0.1)
  # pick most significant proximal and distal clusters
  proximal <- head( filter(sigCluster, label == "proximal") %>% arrange(FDR), 1)
  distal <- head( filter(sigCluster, label == "distal") %>% arrange(FDR), 1)
  status <- c("proximal", "unchanged", "distal", "unchanged")
  status <- data.frame(
    gene = sigGenes[i],
    proximal = "unchanged",
    distal = "unchanged" )

  if( nrow(proximal) > 0){
    status$proximal <- ifelse(proximal$log2FoldChange > 0, "up", "down" )
  }
  if( nrow(distal) > 0){
    status$distal <- ifelse(distal$log2FoldChange > 0, "up", "down" )
  }
  #status <- paste(status, collapse = " ")
  # of the significant clusters, which go up and down?
  return(status)
})

results <- do.call(rbind, clusterDirections)

outFile <- paste0(output.dir, "/", code, "_", mode, "_perGeneDirections.txt")
write.table(results, outFile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  #print(paste( n ,nSig) )
   # if( nrow(cluster) == 0 ){
   #   print(sigGenes[i])
   # }

