library(JunctionSeq) 
library(optparse)
options(echo=T)

option_list <- list(
    make_option(c('--support.tab'), help='', default = ""),
    make_option(c('--gff'), help='', default="/SAN/vyplab/IoN_RNAseq/Kitty/Reference/gencode.vM11.chr_patch_hapl_scaff.basic.annotation.juncseq.gff"),
    make_option(c('--basedir'), help='', default=""),
    make_option(c('--outputdir'), help='', default="")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

support.tab <- opt$support.tab
gff <- opt$gff
basedir <- opt$basedir
outputdir <- opt$outputdir 

#basedir= "/SAN/vyplab/IoN_RNAseq/Kitty/F210I/qorts/" 
#gff="/SAN/vyplab/IoN_RNAseq/Kitty/Reference/gencode.vM11.chr_patch_hapl_scaff.basic.annotation.juncseq.gff" 
#gff="/SAN/vyplab/IoN_RNAseq/Kitty/F210I/qorts/withNovel.forJunctionSeq.gff.gz"
#support.tab = "/cluster/project8/vyp/kitty/F210I/F210I_embryo_June_2016_support.tab"
support = read.table(support.tab,  header = T) 

sample.files <- c() 
for(i in 1:nrow(support)) { 
    sample.files[i] <- paste0(basedir, "qorts/", support$sample[i], "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz") 
} 

print(sample.files) 

# Assumes the sample files have the format 
# sample f1 f2 condition_[cond1] condition_[cond2] 

nconds = ncol(support) - 3 
message("There are ", nconds, " conditions") 

for(i in 1:nconds) { 

   cond.col = i + 3 
   samplestouse = which(!is.na(support[, cond.col])) 
   condition=names(support)[cond.col]

   message("Processing condition: ", condition) 

   if(FALSE){ 
   jscs <- runJunctionSeqAnalyses(sample.files = sample.files[samplestouse],
                sample.names = support$sample[samplestouse],
                condition=support[samplestouse,cond.col],
		flat.gff.file = gff, 
                analysis.type = "junctionsAndExons", nCores=6)

   dir.create(file.path(outputdir, condition)) 
   outfile.prefix=paste0(outputdir, "/", condition, "/juncseq")
   writeCompleteResults(jscs, outfile.prefix=outfile.prefix, save.jscs=TRUE) 

   buildAllPlots(jscs=jscs,
   outfile.prefix = outfile.prefix, 
   use.plotting.device = "svg", 
   FDR.threshold = 0.0001
   )
   } 
   outfile.prefix=paste0(outputdir, "/", condition, "/juncseq")
# Now extract the exons which are significantly DE for downstream analysis 
   jcs.res  <- read.table(paste0(outfile.prefix, "sigGenes.results.txt.gz"), header = T,
                stringsAsFactor = F)


   sig.res <- subset(jcs.res, padjust < 0.01 & featureType == "exonic_part")

   bed.data <- sig.res[, c("chr", "start", "end", "geneID", "padjust", "strand")]

   print(head(bed.data)) 

   out.bed.file <- paste0(outfile.prefix, "sigGenes.bed") 
   out.utr.file <- paste0(outfile.prefix, "sigGenes_utr.bed") 
   write.table(bed.data, quote = F, row.names = F, col.names=F, file = out.bed.file, sep = "\t")

   utr.ref.file = "/SAN/vyplab/IoN_RNAseq/Kitty/Reference/ensembl_3UTR.bed" 
   command=paste0("intersectBed -wa -s -a ", out.bed.file, " -b ", utr.ref.file, " > ", out.utr.file) 
   
   print(command) 
   system(command) 

}  
