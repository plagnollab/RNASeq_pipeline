# junctionseq

Running junctionseq is a two step process. 

Step 1: 
  Run qorts on individual bam files to extract splice junction information 

Step 2: 
  Run the R program junctionseq on the results from qorts 

Requirements:

* BAM files 

* A support file with these columns:
        ``` sample_name bam_file conditionX ```

* A submission file (see examples)

Software requirements: 

* qorts java program 

* junctionseq R package which requires the dexseq package 

 
