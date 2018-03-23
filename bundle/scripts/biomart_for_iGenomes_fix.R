library(biomaRt)


#species <- 'mmusculus_gene_ensembl'
#species <- 'hsapiens_gene_ensembl'
#species <- 'cfamiliaris_gene_ensembl'
#species <- 'drerio_gene_ensembl'
species <- 'sscrofa_gene_ensembl'


if (species == 'xxx') output.file <- paste('zebrafish/biomart/biomart_refseq_', species, '.tab', sep = '')
if (species == 'sscrofa_gene_ensembl') output.file <- paste('pig/biomart/biomart_refseq_', species, '.tab', sep = '')
ensembl <- useMart("ensembl")
ensembl <- useDataset(species,mart=ensembl)
  
attributes <- listAttributes(ensembl)

my.db <- getBM(attributes= c("refseq_mrna", 'ensembl_gene_id', 'wikigene_name', 'chromosome_name', 'strand'),
               mart=ensembl)
my.db <- subset(my.db, refseq_mrna != '')
my.db$strand <- ifelse (my.db$strand == 1, '+', '-')

write.table(x = my.db,
            row.names = FALSE,
            quote = FALSE,
            col.names = FALSE,
            sep = '\t',
            file = output.file)

my.db <- getBM(attributes= c("refseq_mrna_predicted", 'ensembl_gene_id', 'wikigene_name', 'chromosome_name', 'strand'),
               mart=ensembl)
my.db <- subset(my.db, refseq_mrna_predicted != '')
my.db$strand <- ifelse (my.db$strand == 1, '+', '-')

write.table(x = my.db,
            row.names = FALSE,
            quote = FALSE,
            col.names = FALSE,
            sep = '\t',
            file = output.file,
            append = TRUE)




