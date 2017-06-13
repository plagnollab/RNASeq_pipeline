from pybedtools import BedTool
#import csv
import os
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("clusterFile", help = "A bed file of introns that each include an exon")
parser.add_argument("reference", help = "A bed file of iCLIP clusters")
parser.add_argument("outFile", help = "The output file")
args = parser.parse_args()

clusterFile = "/SAN/vyplab/IoN_RNAseq/Nicol/ko/test.bed"
reference = "/SAN/vyplab/HuRNASeq/reference_datasets/mm10.fa"
outFile = "/SAN/vyplab/IoN_RNAseq/Nicol/ko/clean_cluster.bed"

clusterFile = args.clusterFile
reference = args.reference
outFile = args.outFile

for file in [ clusterFile, reference, outFile ]:
	if not os.path.exists(file):
		print str(file) + " does not exist!"
		exit


# for each cluster in clusters grab the coordinates [end+1,end+10] in a strand specific way
# but then due to the clusters being on the opposite strand to the gene, flip the strand. 
# This is necessary to get fasta sequence at the right strand


clusters = BedTool(clusterFile)


def get_downstream_sequence(feature):
	if feature.strand == "-":
		feature.start = feature.end + 1
		feature.end = feature.end + 11
		# switch strand
		feature.strand = "+"
		return(feature)
	elif feature.strand == "+":
		feature.end = feature.start 
		feature.start = feature.start - 10
		# switch strand
		feature.strand = "-"
		return(feature)
	else:
		print "feature must be + or - "
	return(feature)

print("getting flanking coordinates")

downstream = clusters.each(get_downstream_sequence)


# get fasta sequence for each coordinate

print("retrieving genomic sequences")
fasta = downstream.sequence(fi = reference, s = True, tab = True)

# apply tests. 1. is there a stretch of 6*A or are there 7 A in any order?
print("testing each cluster for A-rich sequence")

# random_primers = []
# test = []
# with open(fasta.seqfn, "r") as f:
# 	for line in f:
# 		seq = line.strip().upper().split("\t")[1]
# 		if "AAAAAA" in seq or seq.count("A") >= 7:
# 			random_primers.append(False)
# 		else:
# 			random_primers.append(True)
# 			#print(seq)
# 			test.append(seq)


# # now I have a logical vector of whether the sequence passes or not
# print("writing output")

# #clean_clusters = []
# with open(outFile, "w") as out:
# 	for i in range(len(clusters) ):
# 		if random_primers[i]:
# 			#clean_clusters.append( clusters[i] )
# 			out.write("\t".join(clusters[i]) + "\n")


# do everything in one loop
#outFile2="/SAN/vyplab/IoN_RNAseq/Nicol/ko/clean_cluster_new_method.bed"
# out = open(outFile, "w")
# i = 0
# with open(fasta.seqfn, "r") as f:
# 		for line in f:
# 			if i % 1000 == 0:
# 				print str(i) + " lines processed" 
# 			seq = line.strip().upper().split("\t")[1]
# 			if "AAAAAA" in seq or seq.count("A") >= 7:
# 				i+=1
# 				next
# 			else:
# 				out.write("\t".join(clusters[i]) + "\n")
# 				i+=1
# out.close()


# third method that loops through both the fasta file and the cluster Bed file at the same time:
# much faster
from itertools import izip
out = open(outFile, "w")
i = 0
with open(fasta.seqfn, "r") as sequences, open(clusterFile, "r") as clusters:
	for fasta, clu in izip(sequences,clusters):
		if i % 1000 == 0:
				print str(i) + " lines processed" 
		seq = fasta.strip().upper().split("\t")[1]
		if seq == "":
			print "empty sequence at this position!"
			break 
		if "AAAAAA" in seq or seq.count("A") >= 7:
			next
		else:
			out.write(clu)
		i+=1
out.close()



