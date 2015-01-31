
index=~/Software/misopy-0.4.9/misopy/index_gff.py

#wget http://genes.mit.edu/burgelab/miso/annotations/ver2/miso_annotations_mm10_v2.zip
#unzip miso_annotations_mm10_v2.zip
#mv mm10/ mouse/miso_mm10
#python ~/Software/misopy-0.4.9/misopy/index_gff.py --index mouse/miso_mm10/SE.mm10.gff3  mouse/miso_mm10/indexed_SE_events
python ~/Software/misopy-0.4.9/misopy/index_gff.py --index mouse/miso_mm10/A5SS.mm10.gff3  mouse/miso_mm10/indexed_A5SS_events

#wget http://genes.mit.edu/burgelab/miso/annotations/ver2/miso_annotations_hg19_v2.zip
#unzip miso_annotations_hg19_v2.zip
#mv hg19 human/miso_hg19
