library(data.table)
library(stringr)
library(ggplot2)
library(Biostrings)
# take the output of createVariantTable and make a bunch of BED files

# Fly C9 GR 
sigVarTable <- "/Users/Jack/SAN/HuRNASeq/Fly_C9/SGSeq/control_GR/Fly_C9_splice_variant_table_sig.tab"
nullVarTable <- "/Users/Jack/SAN/HuRNASeq/Fly_C9/SGSeq/control_GR/Fly_C9_splice_variant_table_null.tab"
species <- "Fly" # to specify genome for getFASTA
sgf_object <- "/Users/Jack/SAN/HuRNASeq/Fly_C9/SGSeq/Fly_C9_sgv_novel.RData"

# Nicol FUS KO HOM
sigVarTable <- "/Users/Jack/SAN/IoN_RNAseq/Nicol_FUS/mRNAseq/KO/SGSeq/Control_KO/Nicol_FUS_KO_splice_variant_table_sig.tab"
nullVarTable <- "/Users/Jack/SAN/IoN_RNAseq/Nicol_FUS/mRNAseq/KO/SGSeq/Control_KO/Nicol_FUS_KO_splice_variant_table_null.tab"
species <- "Mouse"
sgf_object <- "/Users/Jack/SAN/IoN_RNAseq/Nicol_FUS/mRNAseq/KO/SGSeq/Nicol_FUS_KO_sgv_novel.RData"



outFolder <- paste0(dirname(sigVarTable), "/bed_files")
if( !dir.exists(outFolder)){dir.create(outFolder)}


createBed <- function(vtable){
  coords <- do.call(rbind, strsplit(vtable$coords,split = ":|-"))
  outBed <- data.frame(
    chr = coords[,1],
    start = coords[,2],
    end = coords[,3],
    gene_name = vtable$gene,
    groupID = vtable$groupID,
    strand = vtable$strand
  )
  return(outBed)
}

if( species == "Fly"){
  fasta <- "/Users/Jack/SAN/HuRNASeq/reference_datasets/dm6.fa"
}
if( species == "Mouse" ){
  fasta <- "/Users/Jack/SAN/HuRNASeq/reference_datasets/mm10.fa"
}
if( species == "Human" ){
  fasta <- "/Users/Jack/SAN/HuRNASeq/reference_datasets/hg38.fa"
}

createBedFASTA <- function(varTable, mode){
  # create bed files for each category of variant
  # mode is either "sig" or "null"
  varTypes <- unique(varTable$variant_type)
  fasta_list <- list()
  bed_list <- list()
  for( v in varTypes ){
    bed <- createBed( varTable[ varTable$variant_type == v,])
    bed_list[[v]] <- bed
    # write out
    bed_out <- paste0(outFolder, "/",v,"_", mode, ".bed")
    fasta_out <- paste0(outFolder, "/",v,"_", mode, ".fasta")
    #print(bed_out)
    write.table(bed, bed_out, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
    # use bedtools to generate fasta sequence
    getFasta.cmd <- paste("bedtools getfasta",
                          "-fi", fasta, 
                          "-bed", bed_out,
                          "-fo", fasta_out,
                          "-s")
    
    system(getFasta.cmd)
    fasta_list[[v]] <- readLines(fasta_out)
  }
  
  # cassette exon central exons 
  central_exons <- paste0(outFolder, "/cassette_exons_central_", mode,".bed")
  if( file.exists( central_exons) ){
    central_exons_df <- read.table(central_exons, header=FALSE, sep = "\t")
    names(central_exons_df) <- c("chr", "start", "end", "geneName", "eventID", "strand")
    bed_list[["central_exons"]] <- central_exons_df
    fasta_out <- gsub(".bed$", ".fasta", central_exons)
    getFasta.cmd <- paste("bedtools getfasta",
                          "-fi", fasta, 
                          "-bed", central_exons,
                          "-fo", fasta_out,
                          "-s")
    
    system(getFasta.cmd)
    fasta_file <- readLines(fasta_out)
    fasta_list[["central_exons"]] <- fasta_file
  }
  
  
  return(list(fasta = fasta_list, bed = bed_list) )
} 

# significant splice variants
sigVarTable <- read.table(sigVarTable, header=TRUE,stringsAsFactors = FALSE)
sig <- createBedFASTA(sigVarTable, "sig")

nullVarTable <- read.table(nullVarTable, header=TRUE, stringsAsFactors = FALSE)
null <- createBedFASTA(nullVarTable, "null")

save(sig,null, file = paste0(outFolder, "/objects.Rdata"))

# TODO: for cassette exons find the internal exon coordinates




# PREDICT FUNCTIONAL CONSEQUENCES FOR THE VARIANT




#### WORKING WITH THE FASTA SEQUENCES


# As the FASTA files are the entire sequence from splice site to splice site I can check for minor splice site usage.
findSpliceSites <- function(fastaList){
  # split off headers (odd) and sequences (even)
  headers <- fastaList[ seq(1, length(fastaList),2)]
  seqs <- str_to_upper(fastaList[ seq(2, length(fastaList),2)])
  # find strand from headers
  strands <- do.call(rbind,strsplit(headers, "\\(|\\)"))[,2]
  spliceSites <- ifelse( strands == "+", 
    yes = paste0(
      str_sub( seqs, start = 1, end = 2), "/",
      str_sub( seqs, start = str_length(seqs) - 2, end = str_length(seqs) -1)
    ),
    no = paste0(
      str_sub( seqs, start = 2, end = 3), "/",
      str_sub( seqs, start = str_length(seqs) -1, end = str_length(seqs))
    )
  )
  return(spliceSites)
} 

sig$spliceSites <- lapply( sig$fasta, FUN = findSpliceSites)
null$spliceSites <- lapply( null$fasta, FUN = findSpliceSites)
# match the splice sites onto the bed files for manual checking

sig$bed$cassette_exon$spliceSite <- sig$spliceSites$cassette_exon
sig$bed$retained_intron$spliceSite <- sig$spliceSites$retained_intron

null$bed$cassette_exon$spliceSite <- null$spliceSites$cassette_exon
null$bed$retained_intron$spliceSite <- null$spliceSites$retained_intron

for( v in names(sig$bed)){
  sig$bed[[v]]$spliceSite <- sig$spliceSites[[v]]
}
for( v in names(null$bed)){
  null$bed[[v]]$spliceSite <- null$spliceSites[[v]]
}


# U12 (minor) introns can also have AG/GT donor and acceptor sequences
# what about the U12 consensus branch point sequence?
# TTCCTTAA should appear within 40nt from the end of the intron
minorBranchSeq <- "TTCCTTA"
#TTCCTTRAY
minorBranch <- function(fastaList){
  # split off headers (odd) and sequences (even)
  headers <- fastaList[ seq(1, length(fastaList),2)]
  seqs <- str_to_upper(fastaList[ seq(2, length(fastaList),2)])
  # find strand from headers
  branchSeqs <- str_sub( seqs, start = str_length(seqs) - 60)
  counts <- str_count( branchSeqs, minorBranchSeq)
  return(counts)
}

sig$minorBranch <- lapply( sig$fasta, FUN = minorBranch)
null$minorBranch <- lapply( null$fasta, FUN = minorBranch)

for( v in names(sig$bed)){
  sig$bed[[v]]$minorBranch <- sig$minorBranch[[v]]
}
for( v in names(null$bed)){
  null$bed[[v]]$minorBranch <- null$minorBranch[[v]]
}

# how many of the sig and null populations have minor branch points?
lapply(null$minorBranch, FUN = function(x){
  prop <- sum(x) / length(x)
  result <- (paste(sum(x), length(x), prop ) ) })

lapply(sig$minorBranch, FUN = function(x){
  prop <- sum(x) / length(x)
  result <- (paste(sum(x), length(x), prop ) ) })


# dinucleotide enrichment 
library(Biostrings)
# waaaay faster
kmerCounts <- function(sequences, k){
  # remove headers
  sequences <- sequences[seq(2, length(sequences), 2)]
  # make one giant string
  sequences <- paste(sequences, collapse = "")
  # use this super optimised function to count
  freq <- oligonucleotideFrequency( DNAString(sequences), k)
  d <- data.frame(freq = freq)
  d$prop <- d$freq / sum(d$freq)
  return(d)
}

kmerCompare <- function( kmer1, kmer2 ){
  #both <- unique( c( kmer1$subseq.total, kmer2$subseq.total))  
  #both <- merge(both, kmer1, by = "subseq.total")
  both <- merge(kmer1, kmer2, by = "row.names")
  #both[ is.na(both)] <- 0
  
    n1 <- sum(both$freq.x)
    n2 <- sum(both$freq.y)
    
    for(i in 1:nrow(both) ){
      both$pvalue[i] <- prop.test(
        x = c(both$freq.x[i], both$freq.y[i]),
        n = c(n1,n2), p = NULL)$p.value
      both$log2FoldChange[i] <- log2(both$prop.x[i]/ both$prop.y[i] )
    }
    both$padj <- p.adjust(both$pvalue, method = "bonferroni")
    names(both)[1] <- "kmer"
    both$kmer <- gsub("T","U",both$kmer)
    both <- both[ order(both$padj),]
  return(both)
}

kmerPlot <- function(kmerComparison, title = NA){
  kmerComparison <- kmerComparison[ kmerComparison$padj < 1,]
  p <- ggplot(kmerComparison, 
              aes( x = log2FoldChange, 
                   y = -log10(padj), 
                   label = kmer )) + 
    geom_text() +
    theme_classic()
  if(!is.na(title)){
    p <- p + labs(title = title)
  }
  return(p)
}



sigRI_6mer <- kmerCounts(sig$fasta$retained_intron, 6)
nullRI_6mer <- kmerCounts(null$fasta$retained_intron, 6)
RI_6mer <- kmerCompare(sigRI_6mer, nullRI_6mer)
kmerPlot(RI_6mer)

sigRI_4mer <- kmerCounts(sig$fasta$retained_intron, 4)
nullRI_4mer <- kmerCounts(null$fasta$retained_intron, 4)
RI_4mer <- kmerCompare(sigRI_4mer, nullRI_4mer)
kmerPlot(RI_4mer)

sigRI_3mer <- kmerCounts(sig$fasta$retained_intron, 3)
nullRI_3mer <- kmerCounts(null$fasta$retained_intron, 3)
RI_3mer <- kmerCompare(sigRI_3mer, nullRI_3mer)
kmerPlot(RI_3mer)

# for each category of splicing variant, compare hexamer enrichment
p_list <- list()
tab_list <- list()
# only first 4 have a reasonable number of events
for( v in names(sig$fasta)[1:4]){
  sigkmer <- kmerCounts(sig$fasta[[v]], 6)
  nullkmer <- kmerCounts(null$fasta[[v]], 6)
  compare <- kmerCompare(sigkmer,nullkmer)
  nEvents <- nrow(sig$bed[[v]])
  title <- paste0(v, " (", nEvents, ")")
  p <- kmerPlot(compare, title = title)
  tab_list[[v]] <- compare
  p_list[[v]] <- p
}
grid.arrange(grobs = p_list, ncol =2)

# merge all together
all_sig <- do.call(c, sig$fasta)
all_null <- do.call(c, null$fasta)
sigkmer <- kmerCounts(all_sig, 6)
nullkmer <- kmerCounts(all_null, 6)
compare <- kmerCompare(sigkmer, nullkmer)
# simulate kmer enrichment by sampling from same group of introns with unbalanced sample size
min_lfc <- c()
min_padj <- c()
for (i in 1:100){
  # sample from null retained introns. Should not be any differences
  nullkmer1 <- kmerCounts( sample(null$fasta$retained_intron[seq(2,length(null$fasta$retained_intron), 2)], 200, replace=FALSE), 3)
  nullkmer2 <- kmerCounts( sample(null$fasta$retained_intron[seq(2,length(null$fasta$retained_intron), 2)], 8000, replace=FALSE), 3)
  compare <- kmerCompare(nullkmer1, nullkmer2)
  min_lfc[i] <- compare[1,]$log2FoldChange
  min_padj[i] <- compare[1,]$padj
}
head(compare)
kmerPlot(compare)

# compare lengths - hypothesis: significantly differentially retained introns are much shorter than null
intronLengths <- function( sequences){
  sequences <- sequences[seq(2, length(sequences), 2)]
  lengths <- str_length(sequences)
  return(lengths)
}


