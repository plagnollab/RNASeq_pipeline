library(biomaRt)

#53     dmelanogaster_gene_ensembl       Drosophila melanogaster genes (BDGP5)

#species <- 'mmusculus_gene_ensembl'
species <- 'hsapiens_gene_ensembl'
#species <- 'cfamiliaris_gene_ensembl'
#species <- 'drerio_gene_ensembl'
#species <- 'sscrofa_gene_ensembl'
#species <- "ggallus_gene_ensembl"
#species <- "dmelanogaster_gene_ensembl"
#species <- "rnorvegicus_gene_ensembl"
#species <- "oaries_gene_ensembl"

output.file <- paste('human/biomart/biomart_annotations_', species, '.tab', sep = '')
ensembl <- useMart("ensembl")
ensembl <- useDataset(species,mart=ensembl)
  
attributes <- listAttributes(ensembl)



my.db <- getBM(attributes= c('ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'start_position', 'end_position', 'strand'),
               mart=ensembl)

names(my.db)[1] <- 'EnsemblID'

if (species == 'hsapiens_gene_ensembl') output.file <- 'human_hg38/biomart/biomart_annotations_human.tab'
if (species == 'mmusculus_gene_ensembl') output.file <- 'mouse/biomart/biomart_annotations_mouse.tab'
if (species == 'cfamiliaris_gene_ensembl') output.file <- 'dog/biomart/biomart_annotations_dog.tab'
if (species == 'drerio_gene_ensembl') output.file <- 'zebrafish/biomart/biomart_annotations_zebrafish.tab'
if (species == 'sscrofa_gene_ensembl') output.file <- 'pig/biomart/biomart_annotations_pig.tab'
if (species == 'ggallus_gene_ensembl') output.file <- 'chicken/biomart/biomart_annotations_chicken.tab'
if (species == 'dmelanogaster_gene_ensembl') output.file <- 'drosophila/biomart/biomart_annotations_drosophila.tab'
if (species == 'rnorvegicus_gene_ensembl') output.file <- 'rat/biomart/biomart_annotations_rat.tab'
if (species == 'oaries_gene_ensembl') output.file <- 'sheep/biomart/biomart_annotations_sheep.tab'



if (species == 'hsapiens_gene_ensembl') {
  my.db <- subset(my.db, chromosome_name %in% c(as.character(1:22), "X", "Y"))
  my.db$chromosome_name <- paste0("chr", my.db$chromosome_name)
}


write.table(x = my.db,
            row.names = FALSE,
            quote = FALSE,
            col.names = TRUE,
            sep = '\t',
            file = output.file)
print(output.file)
