#Make pie chart of the event types in the new datasets 

#library(SGSeq) 
library(dplyr) 
library(ggplot2) 
library(data.table)

library(optparse)
options(echo=T)

# option_list <- list(
#   make_option(c('--support.tab'), help='', default = ""),
#   make_option(c('--step'), help='', default = ""), 
#   make_option(c('--code'), help='', default = ""),
#   make_option(c('--output.dir'), help='', default = ""),
#   make_option(c('--annotation'), help='', default = "")
# )
# 
# 
# option.parser <- OptionParser(option_list=option_list)
# opt <- parse_args(option.parser)
# 
# support.tab <- opt$support.tab
# code <- opt$code
# step <- opt$step 
# output.dir <- opt$output.dir
# species <- opt$species
# annotations.tab <- opt$annotation
# 
# 
# # f210i.dir <- "/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i/"
# # m323k.dir <- "/SAN/vyplab/IoN_RNAseq/Kitty/M323K/sgseq/" 
# 
# f210i.dir <- "/Users/Jack/google_drive/TDP_paper/Sgseq/F210I_embryonic_brain/"
# m323k.dir <- "/Users/Jack/google_drive/TDP_paper/Sgseq/M323K_adult_brain/"
# 
# f210i.res <- paste0(f210i.dir, "F210I_embryonic_brain_res_clean_novel.RData") 
# m323k.res <- paste0(m323k.dir, "M323K_adult_brain_res_clean_novel.RData")
# 
# load(f210i.res) 
# 
# f210i.res.clean <- res.clean 
# 
# load(m323k.res) 
# 
# m323k.res.clean <- res.clean 

#f210i.res.clean <- read.table(f210i.res,header=TRUE)
#m323k.res.clean <- read.table(m323k.res,header=TRUE)

makePieChart <- function(sgseqRes, title, FDRlimit, outFolder = NA){
  # filter by FDR < 0.01!
  res.sig <- dplyr::filter(sgseqRes, FDR < FDRlimit) %>% select(one_of(c("groupID", "variantType", "FDR")) ) 
  # for each groupID take one event
  res.sig.by.group <- res.sig %>% group_by(groupID) %>% do(head(.,1))  
  num.events <- nrow(res.sig.by.group)
  print(paste0("total number of events: ", num.events ))
  print(table(res.sig.by.group$variantType))
  # treat each class of event as separate
  res.events <- table(unlist(strsplit(res.sig.by.group$variantType, "+", fixed= TRUE) ) ) 
  print(res.events)
  print( sum(res.events) )
  res.events.plot <- c() 
  res.events.plot["Cassette exons"] <- sum( c( res.events["S2E:I"], res.events["S2E:S"], res.events["SE:I"], res.events["SE:S"] ), na.rm = TRUE )
  res.events.plot["Retained introns"] <- sum( c( res.events["RI:E"] + res.events["RI:R"] ), na.rm = TRUE )
  res.events.plot["Alternative first exon"] <- sum( c(res.events["AS"] + res.events["AFE"] ), na.rm = TRUE )
  res.events.plot["Alternative last exon"] <- sum( c(res.events["AE"] + res.events["ALE"] ), na.rm = TRUE )
  res.events.plot["Alternative 3' splice site"] <- sum( c(res.events["A3SS:P"] + res.events["A3SS:D"] ), na.rm = TRUE)
  res.events.plot["Alternative 5' splice site"] <- sum( c(res.events["A5SS:P"] + res.events["A5SS:D"] ), na.rm = TRUE)
  res.events.plot["Mutually exclusive exons"] <- sum( c(res.events["MXE"]), na.rm = TRUE )
  res.events.plot <- as.data.frame(res.events.plot) 
  names(res.events.plot) <- "varcounts" 
  print( res.events.plot)
  res.events.plot$variantType <- row.names(res.events.plot)
  res.events.plot <- res.events.plot[ order(res.events.plot$varcounts,decreasing=FALSE),]  
  res.events.plot$variantType <- factor(res.events.plot$variantType, levels = rev(res.events.plot$variantType) )
  res.events.plot <- dplyr::mutate(res.events.plot, pos = cumsum(varcounts) - 0.5*varcounts) 
  res.events.plot$prop <- signif( (res.events.plot$varcounts / sum(res.events.plot$varcounts) ) * 100, 3) 
  res.events.plot$prop <- paste0(res.events.plot$prop, "%")
  res.events.plot <- res.events.plot[ res.events.plot$varcounts >0 ,]
  
  pie <- ggplot(res.events.plot, aes(x="", y = varcounts, fill = variantType)) + 
    geom_bar(width = 1, stat="identity")  + 
    geom_text(aes(x=1.6, y = pos, label = prop), size = 3 ) +
    coord_polar(theta="y") +
    scale_fill_brewer("",palette="Dark2", direction = 1) + 
    theme_void() +
    ggtitle(title) +
    annotate("text", x = 1.8, y = sum(res.events.plot$varcounts) / 2, 
             label = paste0("total number of events at adjusted p < ",FDRlimit,": ", num.events ) )
  
  print(pie)
  if( !is.na(outFolder)){
  ggsave(paste0(outFolder,"/", title, "_sgseq_pie_chart.pdf") )
  
  write.table( res.events.plot, paste0(outFolder, "/", title, "_sgseq_variant_type_table.tab"), col.names = TRUE, sep = "\t")
  }
} 
# 
# makePieChart(f210i.res.clean, "RRM2mut", 0.01, f210i.dir)
# makePieChart(m323k.res.clean, "LCDmut", 0.01, m323k.dir) 


 
