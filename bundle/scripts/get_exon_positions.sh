
gfffile=$1


positionFile=`echo $1 | sed -e 's/gff/position/g'`


echo "EnsemblID exonID chromosome exonStart exonEnd" > $positionFile

grep exonic_part $gfffile | while read chrom name part start end junk1 junk2 junk3 transcripts transcripts_id exonic_part_number exonic_part_numberID gene_id gene_id_id; do
    
    
    enumber=`echo $exonic_part_numberID | sed -e 's/;$//g'`
    echo $gene_id_id $enumber $chrom $start $end >> $positionFile
        
done

echo $positionFile