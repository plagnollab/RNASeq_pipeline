library(Rsubread)
options(echo=T)
species <- 'mouse'

if (species == 'chicken') {
 system('samtools view -b /cluster/scratch3/vyp-scratch2/Daudet_RNASeq/processed/ND1/ND1_unique.bam -o dummy_chicken.bam  1:100-1000')
 gtf <- "chicken/GTF/Gallus_gallus.Galgal4.78.gtf"
 test <- featureCounts  (files = 'dummy_chicken.bam', annot.ext = gtf, isGTFAnnotationFile = TRUE,
                         GTF.featureType = "exon", GTF.attrType = "gene_id")
 
 new.file <- gsub(pattern = ".gtf", replacement = "_length_features.tab", gtf) 
 write.table(x = test$annotation, file = new.file, row.names = FALSE, sep = '\t', quote = FALSE)

}

if (species == 'human') {
  system('samtools view -b /cluster/scratch3/vyp-scratch2/Nejentsev_TB_RNASeq/processed/SLX-7934/A013/A013_unique.bam -o dummy_human.bam  1:100-1000')
  
  test <- featureCounts  (files = 'dummy_human.bam', annot.ext ='human/GTF/human_iGenomes_NCBI37_with_ensembl.gtf', isGTFAnnotationFile = TRUE,
                          GTF.featureType = "exon", GTF.attrType = "gene_id")
  
  write.table(x = test$annotation, file = 'human/GTF/human_iGenomes_NCBI37_with_ensembl_length_features.tab', row.names = FALSE, sep = '\t', quote = FALSE)
}


if (species == 'mouse') {
  gtf <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gtf"
  dummy <- paste0( dirname(gtf), "/dummy_mouse.bam" )
  system( paste0("samtools view -b /cluster/scratch3/vyp-scratch2/IoN_RNASeq/Fratta_RNASeq/MEF/Sample_MEF-M323K-HET1/Sample_MEF-M323K-HET1_unique.bam -o ", dummy, "  1:100-1000") )
  
  test <- featureCounts  (files = dummy, annot.ext = gtf, isGTFAnnotationFile = TRUE,
                          GTF.featureType = "exon", GTF.attrType = "gene_id")
  
  new.file <- gsub(pattern = ".gtf", replacement = "_length_features.tab", gtf)
  write.table(x = test$annotation, file = new.file, row.names = FALSE, sep = '\t', quote = FALSE)


}


if (species == 'Tc1_mouse') {
  system('samtools view -b /cluster/scratch3/vyp-scratch2/IoN_RNASeq/Frances/processed/Tc1_509576/Tc1_509576_unique.bam -o dummy_Tc1mouse.bam  1:100-1000')
  
  test <- featureCounts  (files = 'dummy_Tc1mouse.bam', annot.ext ='Tc1_mouse/GTF/genes.gtf', isGTFAnnotationFile = TRUE,
                          GTF.featureType = "exon", GTF.attrType = "gene_id")
  
  write.table(x = test$annotation, file = 'Tc1_mouse/GTF/Tc1_with_ensembl_length_features.tab', row.names = FALSE, sep = '\t', quote = FALSE)


}


if (species == 'drosophila') {
  system('samtools view -b /cluster/project9/uclgenomics/Bolukbasi/output/Bolukbasi-B2/Bolukbasi-B2_unique.bam -o dummy_drosophila.bam  1:100-1000')
  
  test <- featureCounts  (files = 'dummy_drosophila.bam', annot.ext ='drosophila/GTF/Drosophila_melanogaster.BDGP5.75.gtf', isGTFAnnotationFile = TRUE,
                          GTF.featureType = "exon", GTF.attrType = "gene_id")
  
  write.table(x = test$annotation, file = 'drosophila/GTF/Drosophila_melanogaster.BDGP5.75_length_features.tab', row.names = FALSE, sep = '\t', quote = FALSE)


}
