input=SE.mm10.gff3
output=SE.mm10_NCBI.gff3


sed -e 's/^chr//g' $input | sed -e 's/@chr/@/g' | sed -e 's/=chr/=/g'  > $output