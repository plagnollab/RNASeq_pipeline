library(biomaRt)


species <- 'hsapiens_gene_ensembl'
ensembl <- useMart("ensembl")
ensembl <- useDataset(species,mart=ensembl)
  
attributes <- listAttributes(ensembl)



my.db <- getBM(attributes= c('ensembl_gene_id', 'external_gene_id', 'chromosome_name', 'description', 'start_position', 'end_position'),
               mart=ensembl)

names(my.db)[1] <- 'EnsemblID'

if (species == 'hsapiens_gene_ensembl') output.file <- 'human/biomart/biomart_extra_annotations_human.tab'


my.db$description <- gsub(pattern = ' \\[Source.*', replacement = '', as.character(my.db$description))

write.table(x = my.db,
            row.names = FALSE,
            quote = FALSE,
            col.names = TRUE,
            sep = '\t',
            file = output.file)
print(output.file)
