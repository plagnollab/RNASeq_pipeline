library(biomaRt)

#species <- 'mmusculus_gene_ensembl'
#species <- 'hsapiens_gene_ensembl'
#species <- 'cfamiliaris_gene_ensembl'
#species <- 'drerio_gene_ensembl'
species <- 'sscrofa_gene_ensembl'

output.file <- paste('human/biomart/biomart_annotations_', species, '.tab', sep = '')
ensembl <- useMart("ensembl")
ensembl <- useDataset(species,mart=ensembl)
  
attributes <- listAttributes(ensembl)



my.db <- getBM(attributes= c('ensembl_gene_id', 'external_gene_id', 'chromosome_name', 'start_position', 'end_position', 'strand'),
               mart=ensembl)

names(my.db)[1] <- 'EnsemblID'

if (species == 'hsapiens_gene_ensembl') output.file <- 'human/biomart/biomart_annotations_human.tab'
if (species == 'mmusculus_gene_ensembl') output.file <- 'mouse/biomart/biomart_annotations_mouse.tab'
if (species == 'cfamiliaris_gene_ensembl') output.file <- 'dog/biomart/biomart_annotations_dog.tab'
if (species == 'drerio_gene_ensembl') output.file <- 'zebrafish/biomart/biomart_annotations_zebrafish.tab'
if (species == 'sscrofa_gene_ensembl') output.file <- 'pig/biomart/biomart_annotations_pig.tab'


write.table(x = my.db,
            row.names = FALSE,
            quote = FALSE,
            col.names = TRUE,
            sep = '\t',
            file = output.file)
print(output.file)