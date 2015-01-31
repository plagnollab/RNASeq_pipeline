pythonbin=/share/apps/python-2.7.1/bin/python2.7
pyscript=/cluster/project8/vyp/vincent/libraries/R/installed/DEXSeq/python_scripts/dexseq_prepare_annotation.py

#gtffile=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Mus_musculus/NCBI/build37.2/Annotation/Archives/archive-2012-03-09-05-11-29/Genes/genes.gtf
#gtffile=Mus_musculus.NCBIM37.64.gtf
#${pythonbin} $pyscript $gtffile Mus_musculus.NCBIM37.64.gff
gtffile=/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Dictyostelium_discoideum/GTF/Dictyostelium_discoideum.dictybase.01.23.gtf

#gtffile=$1

gffFile=`echo $gtffile | sed -e 's/gtf/gff/g'`

echo $gtffile $gffFile

${pythonbin} $pyscript $gtffile $gffFile