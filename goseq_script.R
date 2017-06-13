# GOseq script
# Jack Humphrey 2016

# to install all the dependencies
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("goseq","GO.db"))
#install.packages("data.table")

# this loads in the packages
# RUN EACH TIME
library(goseq)
library(data.table)
library(GO.db)
library(optparse)
options(echo=T)
########################## read arguments

option_list <- list(
    make_option(c('--species'), help='', default="mouse"),
    make_option(c('--oFolder'), help='', default = ""),
    make_option(c('--mode'), help='', default="DESeq")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

species <- opt$species

oFolder <- opt$oFolder

mode <- opt$mode



resultsFolder <- paste0(oFolder,"/gene_ontology")
if( !dir.exists(resultsFolder) ){ dir.create(resultsFolder) }
#how to get hold of the deseq results files just from the outFolder? it's predictable though right?
# for each folder present, test if there is a file with the deseq_
dirs <- list.dirs(paste0(oFolder))
print(dirs)
# find all the differential expression results for that dataset
if(mode == "DESeq"){
	results <- list.files(dirs, full.names =T)[ grep("differential_expression.tab",list.files(dirs, full.names =T)) ]
}

if(length(results) == 0){
	stop("no results found!")
}

# this function contains all the instructions to run GOSeq. Feel free to adjust the two thresholds
	## DE.Pvalue.threshold is the initial threshold for uncorrected P values for a gene to be entered into GOSeq. This can be adjusted depending on how strict you want to be.
	## GO.FDR.threshold is the multiple testing threshold for GO enrichment analysis. This doesn't need to be changed.
GOSeq <- function(deseq.res, species, DE.Pvalue.threshold = 0.005, GO.FDR.threshold = 0.05, mode = "DESeq"){
	if(species == "mouse"){
		lengths.file <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed_length_features.tab"
		genome <- "mm10"
	}
	if(species == "human_hg38"){
		lengths.file <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/Homo_sapiens.GRCh38.82_fixed_length_features.tab"
		genome <- "hg38"
	}
	# the output folder is taken to be the same folder that contains the DESeq results file
	#outFolder <- resultsFolder
	#if(length(results) > 1){
		outFolder <- paste0(resultsFolder,"/",basename(dirname(deseq.res)) )
	outFolder <- dirname( deseq.res )
	#
	if(!dir.exists(outFolder) ){ 
		dir.create(outFolder, recursive=T) 
	}


	# read in the data as a dataframe. fread() is from the data.table package and speeds things up a bit
	d <- as.data.frame(fread(deseq.res))
	# remove any genes that do not have an EnsemblID and have a baseMean below 1 - this may be too relaxed. 
		## the idea here is to do the enrichment analysis only with genes that are expressed in some way rather than all the genes known.
	assayed.genes <- subset(d, EnsemblID != "NA" & baseMean >= 1)
	print(paste(nrow(assayed.genes),"out of", nrow(d), "genes are expressed with baseMean > 1"))
	# default mode is DESeq so don't worry.
	if(mode == "DESeq"){
		de.genes <- subset(d, pvalue < DE.Pvalue.threshold & EnsemblID != "NA")$EnsemblID
		de.genes <- unlist(strsplit(de.genes, "+", fixed=T))
		de.genes.number <- paste(length(de.genes), "genes differentially expressed with P <", DE.Pvalue.threshold)
		print(de.genes.number)
	}
	if(mode == "DEXSeq"){
		de.genes <- unique(subset(d, pvalue < DE.Pvalue.threshold & EnsemblID != "NA")$EnsemblID)
		de.genes <- unlist(strsplit(de.genes, "+", fixed=T))
		de.genes.number <- paste(length(de.genes), "genes with exons differentially used with FDR <", DE.Pvalue.threshold)
		print(de.genes.number)
		assayed.genes <- assayed.genes[!duplicated(assayed.genes$EnsemblID),]
	}
	# split apart compound entries
	gene.ids <- assayed.genes$EnsemblID
	gene.ids <- unlist(strsplit(gene.ids, "+", fixed=T))
	# create vector of gene IDs that are differentially expressed
	gene.vector <- as.integer(gene.ids%in%de.genes)
	
	names(gene.vector) <- gene.ids

	# read in gene lengths file and match the gene lengths for each gene.
	lengths <- as.data.frame(fread(lengths.file))
	assayed.lengths <- lengths$Length[match(gene.ids, lengths$GeneID)]
	# use the gene lengths to calculate the length bias and plot it
	pdf(paste0(outFolder,"/bias_graph.pdf"))
	pwf <- nullp(gene.vector, genome,"ensGene", bias.data = assayed.lengths)
	dev.off()

	# actually run the GO analysis
	GO.wall <- goseq(pwf,genome,"ensGene", use_genes_without_cat=FALSE)
	# adjust the p values with Benjamini-Hochberg
	enriched.GO <- GO.wall[p.adjust(GO.wall$over_represented_pvalue, method="BH")<GO.FDR.threshold,]
	
	# write out the results plus an explanation of each significant GO term.
	GO.res <- paste0(outFolder,"/enriched_GO_report.txt")
	GO.exp <- paste0(outFolder,"/enriched_GO_terms_explanation.txt")

	# if nothing is found
	if(nrow(enriched.GO) == 0){
		null.statement <- paste("no significant GO terms")
		print(null.statement)
		write.table(null.statement, GO.res, row.names=F,quote=F,sep="\t")
		return()
	}

	# get the GO annotations for each DE gene
	go <- getgo(names(gene.vector[gene.vector == 1]),genome,"ensGene")
	terms <- enriched.GO$category
	terms_human <- paste(enriched.GO$category, enriched.GO$term)
	GO.genes <- list()
	for(i in 1:length(terms)){GO.genes[[i]] = names(gene.vector[gene.vector == 1])[grep(terms[i],go)] }
	names(GO.genes) <- terms_human
	# deal with compound genes
	de.gene.names <- subset(d, pvalue < DE.Pvalue.threshold & EnsemblID != "NA")$external_gene_id
	de.gene.names <- unlist(strsplit(de.gene.names, "+", fixed=T))


	GO.genes.human <- lapply(X = GO.genes, FUN = function(x) de.gene.names[match(x,de.genes)]) 
	d$direction <- ifelse(d$log2FoldChange > 0, "UP","DOWN")
	GO.genes.human <- lapply(X = GO.genes.human, FUN = function(x) paste(x, d$direction[match(x, d$external_gene_id)]))

	GO.genes.out <- paste0(outFolder,"/GO_term_genes_ensembl.txt")
	GO.genes.human.out <- paste0(outFolder,"/GO_term_genes_human.txt")
	
	write(sapply(names(GO.genes),function(x) paste(x,"\t",paste(GO.genes[[x]],collapse="\t"))), file = GO.genes.out )
	write(sapply(names(GO.genes.human),function(x) paste(x,"\t",paste(GO.genes.human[[x]],collapse="\t"))), file = GO.genes.human.out )
 	
	if(nrow(enriched.GO) > 0){
		# write results to a table
		write.table(enriched.GO, GO.res, row.names=F,quote=F,sep="\t")
		
		print(paste(nrow(enriched.GO),"GO terms significantly over-represented with FDR<", GO.FDR.threshold))
		# print out details of all (or first 25) enriched GO terms
		sink(GO.exp)
		print(paste(length(de.genes),"genes differentially expressed with P<", DE.Pvalue.threshold))
		print(paste(nrow(enriched.GO),"GO terms significantly over-represented with FDR<", GO.FDR.threshold))
		for(go in enriched.GO$category[1:length(enriched.GO$category)]){
			print(GO.wall[GO.wall$category == go,]$over_represented_pvalue)
			print(GOTERM[[go]])
			cat("--------------------------------------\n")
		}

		sink()
		
	}
}

for(i in results){
	GOSeq(deseq.res = i, species = species, mode = mode)
}


quit()

# FOR TESTING


DE.Pvalue.threshold = 0.005
GO.FDR.threshold = 0.05
mode = "DESeq"
# run the function with your results
GOSeq(deseq.res,species)
#Change these values for your own analysis
	##deseq.res is your differential expression results file from DESeq2
deseq.res <- "/SAN/vyplab/HuRNASeq/Baloh_C9_iPSC/processed/deseq2/control_MND_c9/deseq_baloh_differential_expression.tab"
	##define your species for the gene lengths file. This only works for human build 38 and mouse build 10.
species <- "human_hg38"

deseq.res <- "/SAN/vyplab/IoN_RNAseq/Anny_FUS/12_months/Spinal_Cord/deseq2/FUS_WT_FUS_D14/deseq_Anny_12months_Spinal_Cord_differential_expression.tab" 
species <- "mouse"

dupuis <- c('/SAN/vyplab/HuRNASeq/Dupuis_FUS_NLS/FUS_NLS/deseq2/CTL_FUS_deltaNLS/deseq_Dupuis_NLS_differential_expression.tab',
		'/SAN/vyplab/HuRNASeq/Dupuis_FUS_NLS/FUS_KO/deseq2/CTL_FUS_KO/deseq_Dupuis_KO_differential_expression.tab')

cleveland_fus <- "/SAN/vyplab/HuRNASeq/Cleveland_FUS/processed/deseq2/CTL_FUS/deseq_Cleveland_FUS_differential_expression.tab"
