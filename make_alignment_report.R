library(tidyr)
library(stringr)

args <- commandArgs(trailingOnly=TRUE)

#print(args)


outFolder <- args[1]

#print(outFolder)

all_logs <- list.files(outFolder, pattern = "Log.final.out", recursive = TRUE, full.names = TRUE )

print(all_logs)

logbook <- list()

for( i in 1:length(all_logs) ){
	log <- readLines(all_logs[i])
	
	lines <- str_trim( log[grepl("niquely",log) ] )
	lines <- str_split_fixed(lines, "\\|\t", 2)
	sample <- basename(dirname(all_logs[i]))
	logbook[[i]] <- data.frame( sample = sample, key = lines[,1], value = lines[,2] )
}

df <- do.call(what = rbind, args = logbook )

df <- spread(df, "key","value")

# format

names(df) <- c("Sample", "Unique_mapping_%", "Unique_mapping_number")

df$Unique_mapping_number <- format( as.numeric( str_trim(df$Unique_mapping_number) ), big.mark = "," )

write.table( df, paste0(outFolder,"/alignment_report.tab"), sep = "\t", quote = FALSE , col.names = TRUE, row.names = FALSE )

