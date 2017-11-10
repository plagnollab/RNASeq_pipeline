# createVariantTable
# collapse complicated SGSeq output into a table of splice variants
library(dplyr)
library(stringr)
library(ggplot2)
# for testing
sgseq_res <- "/Users/Jack/SAN/IoN_RNAseq/FTD_brain/SGSeq/CTL_FTD_TAU/FTD_brain_CTL_FTD_TAU_res_clean_novel.tab"
support_frame <- "/Users/Jack/SAN/IoN_RNAseq/FTD_brain/SGSeq_support.tab"
condition <- "condition_MAPT"

# FLY + GR dipeptides
sgseq_res <- "/Users/Jack/SAN/HuRNASeq/Fly_C9/SGSeq/control_GR/Fly_C9_control_GR_res_clean_novel.tab"
support_frame <- "/Users/Jack/SAN/HuRNASeq/Fly_C9/fly_SGSeq_support.tab"
code <- "Fly_C9_GR"
condition <- "conditionGR"
# Fly + GA dipeptides
sgseq_res <- "/Users/Jack/SAN/HuRNASeq/Fly_C9/SGSeq/control_GA/Fly_C9_control_GA_res_clean_novel.tab"
condition <- "conditionGA"
code <- "Fly_C9_GA"

# Nicol FUS KO
sgseq_res <- "/Users/Jack/SAN/IoN_RNAseq/Nicol_FUS/mRNAseq/KO/SGSeq/Control_KO/Nicol_FUS_KO_Control_KO_res_clean_novel.tab"
support_frame <- "/Users/Jack/SAN/IoN_RNAseq/Nicol_FUS/mRNAseq/KO/nicol_KO_SGSeq_support.tab"
code <- "Nicol_FUS_KO"
condition <- "condition_HOM"
sgf_object <- "/Users/Jack/SAN/IoN_RNAseq/Nicol_FUS/mRNAseq/KO/SGSeq/Nicol_FUS_KO_sgv_novel.RData"

# Nicol FUS d14
sgseq_res <-"/Users/Jack/SAN/IoN_RNAseq/Nicol_FUS/mRNAseq/d14/SGSeq/Control_HOM/Nicol_FUS_d14_Control_HOM_res_clean_novel.tab"
support_frame <- "/Users/Jack/SAN/IoN_RNAseq/Nicol_FUS/mRNAseq/d14/nicol_d14_SGSeq_support.tab"
code <- "Nicol_FUS_d14"
condition <- "condition_HOM"

# BEGIN

d <- read.table(sgseq_res, header=TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
support <- read.table(support_frame, header=TRUE, sep = "\t", stringsAsFactors = FALSE)
# get samples from support
support.loc <- support[, c("sample_name", condition)]
support.loc <- support.loc[ complete.cases(support.loc),]
sample_list <- list()
for( cond in unique(support.loc[,condition]) ){
  sample_list[[cond]] <- support.loc$sample_name[ support.loc[, condition ] == cond ]
}


variantType_df <- data.frame(
  SGSeq = c("SE", "S2E", "RI", "A3SS", "A5SS", "ALE", "AFE", "AE", "AS", "MXE"),
  human = c("cassette_exon", "multi-cassette_exon", "retained_intron", "alt_3_splice_site", "alt_5_splice_site", "alt_last_exon", "alt_first_exon", "multi-alt_last_exon", "multi_alt_first_exon", "mut_exclusive_exons"),
  stringsAsFactors = FALSE
)

##########
# FUNCTIONS
##########

# create a cleaned up table from the SGSeq results which should be easier to understand
# as SGSeq treats each event as binary but gives each outcome PSI values
# a reference outcome needs to be set and this can be based on the type of event
# with casette exons the outcome of exon inclusion is the reference outcome
# dPSIs are calculated for the reference outcome
createVarTable <- function(d, groupIDs){
  varTable <- lapply(groupIDs, FUN = function(i){
    event <- d[ d$groupID == i,]
    varType <- event$variantType
    varTypeSplit <- str_split_fixed(varType, ":",2)[,1]
    varTypeHuman <- NA
    # simple variants first:
    if( !any(grepl("\\+", varType)) ){
      ref <- "unassigned"
      if( nrow(event) == 1){
        ref <- 1
      }else{
        # for each event, rule which of the two is the reference
        # if any other type (ALE, AFE, MXE etc) or a mixture then just pick the most represented (meanBase)
        if( all("SE" %in% varTypeSplit) ){
          ref <- which( varType == "SE:I")
        }
        if( all("S2E" %in% varTypeSplit) ){
          ref <- which( varType == "S2E:I")
        }
        if( all("RI" %in% varTypeSplit) ){
          ref <- which( varType == "RI:R")
        }
        if( all("A5SS" %in% varTypeSplit) ){
          ref <- which( varType == "A5SS:P")
        }
        if( all("A3SS" %in% varTypeSplit) ){
          ref <- which( varType == "A3SS:P")
        }
        if( length(ref) > 1 | ref == "unassigned"){
          ref <- which( event$exonBaseMean == max(event$exonBaseMean))
        }
      }
    }else{
      # for complex variants 
      #set whichever row has the highest meanBase as ref
      ref <- head(which( event$exonBaseMean == max(event$exonBaseMean)), 1)
      #return( list( i, varType, ref, "complex"))
      varTypeHuman <- "complex"
    }
    # create dPSI from the reference row
    control_psi <- event[ref, paste0(sample_list[[1]], "_psi") ]
    case_psi <- event[ref, paste0(sample_list[[2]], "_psi") ]
    dPSI <- mean(as.numeric(case_psi), na.rm = TRUE) - mean(as.numeric(control_psi), na.rm = TRUE)
    
    # translate varType to understandable
    if(is.na(varTypeHuman) ){
      if( !any(grepl("\\+", varType)) ){
        varTypeHuman <- variantType_df$human[ match( unique(varTypeSplit), variantType_df$SGSeq)]
        if( all(is.na(varTypeHuman) ) | length(varTypeHuman) > 1 ){
          varTypeHuman <- "complex"
        }
      }else{
        varTypeHuman <- "complex"
      }
    }
    # pick the largest coordinates to encompass all variants
    coords <- do.call(rbind,strsplit(event$coords, split = ":|-"))
    coordSizes <- apply(coords, 1, FUN = function(x) as.numeric(x[3]) - as.numeric(x[2]) )
    biggestCoords <- head(event$coords[order(coordSizes)],1)
    
    # annotation - are both the reference and alternate forms annotated?
    ref_anno <- event$txName[ref]
    alt_anno <- event$txName[-ref]
    if( all(ref_anno == "") ){ ref_anno <- "novel"}else{ ref_anno <- "annotated"}
    if( all(alt_anno == "") ){ alt_anno <- "novel"}else{ alt_anno <- "annotated"}
    #ref_anno <- paste(ref_anno, collapse = "+")
    #alt_anno <- paste(alt_anno, collapse = "+")
    # if dPSI > 0 then the ref form is increasing - is it annotated?
    # if dPSI < 0 then the alt form is increasing - is it annotated?
    # if( dPSI > 0){
    #   annotation <- refAnno
    # }else{
    #   annotation <- altAnno
    # }
    # annotation <- ifelse(annotation == "", yes = "novel", no = annotation)
    # stranding! This should have been created in step2b oh well
    strand <- unique( str_split_fixed(c(event$to,event$from), ":", 4)[,4] )
    if( length(strand) > 1){
      strand <- "*"
    }
    # assemble into a table
    varTable_row <- data.frame(
      groupID = i,
      gene = event$geneName[ref], # in case of weirdness
      EnsemblID = event$ensemblName[ref],
      coords = biggestCoords,
      strand = strand,
      variant_type = varTypeHuman,
      control_PSI = signif( mean(as.numeric(control_psi), na.rm = TRUE), 3),
      case_PSI = signif( mean(as.numeric(case_psi), na.rm = TRUE), 3),
      dPSI = signif(dPSI,3),
      FDR = signif(min(event$padj, na.rm= TRUE),3),
      ref_anno,
      alt_anno,
      stringsAsFactors=FALSE
    )
    return(varTable_row)
  })
  varTable <- do.call(rbind,varTable)
return(varTable)
}

# plot distribution of delta PSI by variant type
# to do: create a Venn
dPSIPlot <- function(varTable, outFile){
  variants <- table(varTable$variant_type)
  variants <- names(variants[ order(variants, decreasing = TRUE)])
  varTable$variant_type <- factor( varTable$variant_type, levels = variants)
  
  p <- ggplot(varTable, 
              aes(y = dPSI, 
                  x = variant_type, 
                  group = variant_type
              )) + 
    geom_jitter(aes(size = -log10(FDR), colour = variant_type)) +
    geom_boxplot(fill = NA, outlier.colour = NA, width = 0.5 ) +
    theme_classic() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    guides(colour=FALSE) +
    theme(axis.text.x  = element_text(angle=45, hjust=1) ) +
    xlab("")
  print(p)
  ggsave(p,filename = gsub(outFile,".tab", ".pdf") )
}


# get all significant group IDs - sometimes only one of the events will be significant!
sigOutFile <- paste0(dirname(sgseq_res), "/", code ,"_splice_variant_table_sig.tab")
sigGroupIDs <- unique( dplyr::filter(d, padj < 0.05 & geneName != "")$groupID)
print("Number of unique events:")
print(length(sigGroupIDs))

sigVarTable <- createVarTable(d, sigGroupIDs)
write.table(sigVarTable, sigOutFile, sep = "\t", row.names = FALSE, quote = FALSE)

dPSIPlot(sigVarTable, sigOutFile)

# do the same but for non-significant splicing events
nullOutFile <- paste0(dirname(sgseq_res), "/", code ,"_splice_variant_table_null.tab")
nullGroupIDs <- unique( dplyr::filter(d, padj > 0.95 & geneName != "")$groupID)
print("Number of unique null events:")
print(length(nullGroupIDs))
# create table
nullVarTable <- createVarTable(d, nullGroupIDs)
write.table(nullVarTable, nullOutFile, sep = "\t", row.names = FALSE, quote = FALSE)

dPSIPlot(nullVarTable, nullOutFile)


# from Kitty

load(sgf_object)

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




# predict effect of cassette splicing variants
if( species == "mouse"){
  genome <- "BSgenome.Mmusculus.UCSC.mm10"
  transcripts <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
}

library(genome, character.only = TRUE)
library(transcripts, character.only = TRUE)

txdb <- eval(parse(text = transcripts))
genome <- eval(parse(text = genome))
load(sgf_object)

test <- sgv_novel[ which(SGSeq::eventID(sgv_novel) == 38578 ) ]
vep <- predictVariantEffects(test, tx = txdb, genome = genome, output = "full", cores = 2)

# vep is a dataframe - can be manipulated
