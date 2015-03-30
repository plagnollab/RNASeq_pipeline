library(Rsubread)

species <- 'chicken'

if (species == 'chicken') {
 system('samtools view -b /scratch2/vyp-scratch2/Daudet_RNASeq/processed/ND1/ND1_unique.bam -o dummy_chicken.bam  1:100-1000')
 gtf <- "chicken/GTF/Gallus_gallus.Galgal4.78.gtf"
 test <- featureCounts  (files = 'dummy_chicken.bam', annot.ext = gtf, isGTFAnnotationFile = TRUE,
                         GTF.featureType = "exon", GTF.attrType = "gene_id")
 
 new.file <- gsub(pattern = ".gtf", replacement = "_length_features.tab", gtf) 
 write.table(x = test$annotation, file = new.file, row.names = FALSE, sep = '\t', quote = FALSE)

}

if (species == 'human') {
  system('samtools view -b /scratch2/vyp-scratch2/Nejentsev_TB_RNASeq/processed/SLX-7934/A013/A013_unique.bam -o dummy_human.bam  1:100-1000')
  
  test <- featureCounts  (files = 'dummy_human.bam', annot.ext ='human/GTF/human_iGenomes_NCBI37_with_ensembl.gtf', isGTFAnnotationFile = TRUE,
                          GTF.featureType = "exon", GTF.attrType = "gene_id")
  
  write.table(x = test$annotation, file = 'human/GTF/human_iGenomes_NCBI37_with_ensembl_length_features.tab', row.names = FALSE, sep = '\t', quote = FALSE)
}


if (species == 'mouse') {
  system('samtools view -b /scratch2/vyp-scratch2/Fratta_RNASeq/MEF/Sample_MEF-M323K-HET1/Sample_MEF-M323K-HET1_unique.bam -o dummy_mouse.bam  1:100-1000')
  
  test <- featureCounts  (files = 'dummy_mouse.bam', annot.ext ='mouse/GTF/mouse_iGenomes_GRCm38_with_ensembl.gtf', isGTFAnnotationFile = TRUE,
                          GTF.featureType = "exon", GTF.attrType = "gene_id")
  
  write.table(x = test$annotation, file = 'mouse/GTF/mouse_iGenomes_GRCm38_with_ensembl_length_features.tab', row.names = FALSE, sep = '\t', quote = FALSE)


}


if (species == 'Tc1_mouse') {
  system('samtools view -b /scratch2/vyp-scratch2/IoN_RNASeq/Frances/processed/Tc1_509576/Tc1_509576_unique.bam -o dummy_Tc1mouse.bam  1:100-1000')
  
  test <- featureCounts  (files = 'dummy_Tc1mouse.bam', annot.ext ='Tc1_mouse/GTF/genes.gtf', isGTFAnnotationFile = TRUE,
                          GTF.featureType = "exon", GTF.attrType = "gene_id")
  
  write.table(x = test$annotation, file = 'Tc1_mouse/GTF/Tc1_with_ensembl_length_features.tab', row.names = FALSE, sep = '\t', quote = FALSE)


}
