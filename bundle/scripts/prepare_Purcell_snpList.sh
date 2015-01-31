#awk '{print $3}' human/exonic_snps_human_hg19.tab | sed -e 's/:/ /g' > human/exonic_snps_human_hg19_clean.tab

awk '{print $3}' human/exonic_snps_human_hg19.tab | sed -e 's/:/ /g' | awk '{print $1":"$2-2"-"$2+2}' > human/exonic_snps_human_hg19_regions.tab