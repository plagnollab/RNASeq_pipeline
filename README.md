Set of scripts for RNA-Seq data processing, in particular differential expression analysis

# Description of pipeline

* After optional adapter trimming with Trim Galore! (0.4.1), the fastq reads were aligned using STAR (2.5.0a_alpha) to the human build 38 (GCA_000001405.15_GRCh38_no_alt_analysis_set).
* The resulting BAM file was sorted using NovoSort(1.00.01) and duplicate reads were flagged using Picard MarkDuplicates(1.100).
* The GRCh38 GTF transcript file was flattened to create a GFF file, a set of union exons, using the dexseq_prepare_annotation.py Python script included with the DEXSeq package.
* The aligned reads overlapping the union exons were counted using HTSeq. This is wrapped in the dexseq_count.py Python script, included with the DEXSeq package.
* Differential exon and transcript expression between conditions was assessed using the DEXSeq (1.14.2) and DESeq2 (1.8.2) Bioconductor packages respectively running on R (3.1.1).


# Requirements

[![Join the chat at https://gitter.im/plagnollab/RNASeq_pipeline](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/plagnollab/RNASeq_pipeline?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

R packages and software:

- R = 3.2.2
- BioConductor 3.1
```
/cluster/project8/vyp/vincent/Software/R-3.2.2/bin/R
```
- DESeq2  >= 1.8.2
- DEXSeq >= 1.14.2
- BiocParallel >= 1.2.22 
- data.table
- optparse

Notes for installation of DEXSeq/DESeq2:

- Need to add /share/apps/binutils-2.25/bin to the PATH
- type 
```
scl enable devtoolset-1.1 'bash'
```
and then run R

# Input format

The key input is a table in plain text format that contains one row per sample, with the header line:
```
sample f1 f2 condition
```
Assuming that condition is what you want to run the differential expression analysis on.
You also need to specify the input folder, so that the fastq can be found at ${iFolder}/${f1} and ${iFolder}/${f2}.
Note that f1 and f2 can specify subfolders themselves. If the data is single stranded then the f2 column should be present but have the value of NA for each sample. 
Also a species parameter and and output folder.


