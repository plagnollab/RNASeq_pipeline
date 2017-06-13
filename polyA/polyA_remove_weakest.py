# polyA tools
# remove any clusters that contribute < 5% of the total number of reads to its host gene

from pybedtools import BedTool
#import csv
import os
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("clusterFile", help = "A bed file of introns that each include an exon")
parser.add_argument("minProportion", help = "minimum contribution a single cluster can make to the total counts")
parser.add_argument("outFile", help = "The output file")
args = parser.parse_args()

clusterFile = args.clusterFile
minProportion = float(args.minProportion)
outFile = args.outFile

#minProportion = 0.05

#clusterFile = "/SAN/vyplab/IoN_RNAseq/Nicol/ko/test.bed"
#outFile = "/SAN/vyplab/IoN_RNAseq/Nicol/ko/test_thinned.bed"

clusters = BedTool(clusterFile)

# create a dictionary storing the clusters for each gene
genes = {}
for feature in clusters:
	gene = str(feature.fields[6])
	if gene not in genes:
		genes[gene] = [ feature ]
	else:
		genes[gene].append( feature )


# go through and remove any that contribute less than 5% of the total counts
with open(outFile, "w") as out:
	for gene in genes:
		total = 0
		#print gene
		for cluster in genes[gene]:
			#print cluster.score
			#total+=cluster.score
			total += int(cluster.name)
		#print(total)

		# now go back and remove any 
		threshold = int(minProportion * total)
		#print "threshold:", threshold
		#print "length:", len(genes[gene])
		for cluster in genes[gene]:
			#print cluster.name
			if int(cluster.name) < threshold:
				#print "cluster thrown out!"
				genes[gene].remove(cluster)
			else:
				out.write("\t".join(cluster) + "\n")
		#print "new length:", len(genes[gene])

# write out to file

