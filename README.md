Set of scripts for RNA-Seq data processing, in particular differential expression analysis

# Description of pipeline (RNA-Seq pipeline version 8)

* After optional adapter trimming with Trim Galore! (0.4.1), the fastq reads were aligned using STAR (2.4.2a) to the human build 38 (GCA_000001405.15_GRCh38_no_alt_analysis_set).
* The resulting BAM file was sorted and duplicate reads were flagged using NovoSort(1.03.09).
* The GRCh38 GTF transcript file was flattened to create a GFF file, a set of union exons, using the dexseq_prepare_annotation.py Python script included with the DEXSeq package.
* The aligned reads overlapping the union exons were counted using HTSeq. This is wrapped in the dexseq_count.py Python script, included with the DEXSeq package.
* Differential exon and transcript expression between conditions was assessed using the DEXSeq (1.14.2) and DESeq2 (1.8.2) Bioconductor packages respectively running on R (3.1.1).

<p align="center">
  <img src="https://github.com/plagnollab/RNASeq_pipeline/blob/master/schematic.png">
</p>

# Requirements

[![Join the chat at https://gitter.im/plagnollab/RNASeq_pipeline](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/plagnollab/RNASeq_pipeline?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
  
- STAR (>= 2.4.2a)
- NovoSort(>= 1.03.09)
  
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

# Input 1: Sample Table

The key input is a table in **tab delimited** plain text format that contains one row per sample, with the header line:
```
sample  f1  f2  condition  (type)
```

* `sample` must be a unique sample name.
* `f1` and `f2` specify the input fastqs, which can be gzipped. If the data are single stranded then f2 must be set as **NA**.
* If the fastq files are in pieces then fastq files of the same direction should be comma separated.
* `condition` is what you want to run the differential expression analysis on. You can include multiple condition columns
* `type` is an optional covariate column which will be included in the differential expression model.
* The input folder is specified in the submission form. However if the fastq files are in separate subfolders then these subfolders should be present in the sample table, ie `/subfolder/sample1.fastq`.

# Input 2: Submission Form
The second input file is a list of variables that will be used by the pipeline script. An example is included in the repository and should be filled in by the user. Each variable is named and explained below:
* `pipeline`: path to where the repository is downloaded
* `oFolder`: where the processed data will go. Conventionally `samples/processed`.
* `iFolder`: where the input fastq files are. If in separate subfolders then include the subfolder in the sample table.
* `dataframe`: the path to the sample table
* `code`: the name of the project, used as the job name when submitting each step.
* `species`: the species genome used for alignment. `human` specifies hg19 whereas `human_hg38` specifies hg38.
* `stranded=(unstranded/fr-firststrand/fr-secondstrand)`: whether the RNA-seq library is stranded. 
* `submit=(yes|no)`: whether the pipeline should automatically submit jobs to the cluster.
* `step0_QC=(yes|no)`: each fastqc is checked with FastQC. This step should ideally be completed before the rest of the pipeline is run.
* `trim_galore=(yes|no)`: should the fastq files be trimmed before alignment?
* `starStep1=(yes|no)`: the initial alignment step with STAR, followed by sorting and duplicate marking.
* `starStep1b=(yes|no)`: indexing of the resulting bam file and computing of summary statistics.
* `starStep2=(yes|no)`: counting the reads in each sample with HTSeq.
* `prepareCounts=(yes|no)`: preparation of read counts for DESeq and DEXSeq; creation of per-sample FPKMs.
* `Rdeseq=(yes|no)`: differential gene expression with the DESeq2 package.
* `Rdexseq=(yes|no)`: differential exon usage with the DEXSeq package.

## A note on stranding
Strand-specific library preparations are now commonplace, which improves the accuracy of feature quantification. However, there are two possible ways of stranding a paired RNAseq library:
* fr-firststrand, where the first read of the pair maps to the same orientation of the gene.
* fr-secondstrand, where the first read of the pair maps to the opposite orientation of the gene.
If you're unsure whether the library is forward or reverse stranded, set stranded=unstranded and run step 1a to align your reads. Then run infer_experiment.py from the RSeQC package (http://rseqc.sourceforge.net/#infer-experiment-py).

## Advanced usage
Two flags in the submission script, `summary` and `force` are now deprecated. They have now been repurposed for non-standard use cases.
#### Only outputting lists of splice junctions from STAR
* Set `starStep1a` to `yes` and `force` to `SJsOnly`.

#### Realign files that have already been trimmed
  Useful if pipeline has crashed downstream of trimming.
* Set `trim_galore` to `yes` and `summary` to `trimmed_exist`.

#### Two-pass mapping with STAR
  If you're interested in novel/unannotated splicing events you should consider using STAR's two-pass mapping mode. 
  This first aligns your reads using a known set of transcripts (provided by Ensembl), with any novel splice junctions being mapped. It then incorporates those novel junctions into the referene and re-aligns.
  The two-pass approach will not necessarily increase the detection of novel junctions but will improve the number of splice reads mapping to them.
  It is however much slower than the usual mode.
  * Set `summary` to `twopass`

#### Chimeric alignments
  If you're interested in circular RNAs or fusion transcripts then you can turn on STAR's reporting of chimeric alignments with
* Set `force` to `chimeric`
  This will align samples as usual but for each sample will output separate `Chimeric.out.sam` and `Chimeric.out.junction` files. See the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) for more details.
