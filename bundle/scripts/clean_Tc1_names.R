options(stringsAsFactors = FALSE)
data <- read.table('Tc1_mouse/GTF/Tc1_length_features_old.tab', header = TRUE)
mouse <- read.table('mouse/biomart/biomart_annotations_mouse.tab', header = TRUE)
human <- read.table('human/biomart/biomart_annotations_human.tab', header = TRUE)


new.name <- data$GeneID
new.name <- ifelse ( new.name %in% mouse$external_gene_id, mouse$EnsemblID [ match(new.name, mouse$external_gene_id) ], new.name)
new.name <- ifelse ( new.name %in% human$external_gene_id, human$EnsemblID [ match(new.name, human$external_gene_id) ], new.name)

data$GeneID <- new.name

write.table(x = data, file = 'Tc1_mouse/GTF/Tc1_length_features.tab', row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
