#!/share/apps/perl-5.14.2/bin/perl

#use strict;   

my $species = $ARGV[0];

my $biomart;
my $output;
my $inputGTF;
my $transcriptNB = 11;

if ($species eq "dog") {
    $biomart = "dog/biomart/biomart_refseq_cfamiliaris_gene_ensembl.tab";
    $output = "dog/GTF/dog_iGenomes_NCBI_3_1_with_ensembl.gtf";
    $inputGTF = "/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Canis_familiaris/NCBI/build3.1/Annotation/Genes/genes.gtf";
    $transcriptNB = 15
}

if ($species eq "pig") {
    $biomart = "pig/biomart/biomart_refseq_sscrofa_gene_ensembl.tab";
    $output = "pig/GTF/pig_iGenomes_NCBI_10_2_with_ensembl.gtf";
    $inputGTF = "/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Sus_scrofa/NCBI/Sscrofa10.2/Annotation/Genes/genes.gtf";
    $transcriptNB = 15
}

if ($species eq "mouse") {
    $biomart = "mouse/biomart/biomart_refseq_mmusculus_gene_ensembl.tab";
    $output = "mouse/GTF/mouse_iGenomes_GRCm38_with_ensembl.gtf";
    $inputGTF = "/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Mus_musculus/NCBI/GRCm38/Annotation/Genes/genes.gtf";
    $transcriptNB = 15
}


if ($species eq "zebrafish") {
    $biomart = "zebrafish/biomart/biomart_refseq_drerio_gene_ensembl.tab";
    $output = "zebrafish/GTF/zebrafish_iGenomes_Zv9_with_ensembl.gtf";
    $inputGTF = "/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Danio_rerio/NCBI/Zv9/Annotation/Genes/genes.gtf";
    $transcriptNB = 15
}



if ($species eq "human") {
    $biomart = "human/biomart/biomart_refseq_hsapiens_gene_ensembl.tab";
    $output = "human/GTF/human_iGenomes_NCBI37_with_ensembl.gtf";
    $inputGTF = "/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Homo_sapiens/NCBI/build37.2/Annotation/Genes/genes.gtf";
}

if ($species eq "human_muscle") {
    $biomart = "support/biomart_refseq_hsapiens_gene_ensembl.tab";
    $output = "humanmuscle_iGenomes_NCBI37_with_ensembl.gtf";
    $inputGTF = "/SAN/biomed/biomed14/vyp-scratch/vincent/tophat_reference/Homo_sapiens/NCBI/build37.2/Annotation/Genes/genes.gtf";
}

print "Looking at species $species\n";

print "Reading the file $biomart\n";
open (INP1, " < $biomart") or die "Cannot read $biomart";
my %hash = ();
my %geneToChrom = ();
my %transcriptToChrom = ();
my %geneToStrand = ();


print "Output file in $output\n";

while (<INP1>) {
    chomp $_;
    my @spl1 = split('\t', $_);


    ##In human muscle, if the gene is DMD, only keep NM_004006
    if (($species eq "human_muscle") && ($spl1[ 2 ] eq "DMD")) {
	if ($spl1[0] ne "NM_004006") {next;}
    }

    $hash{ $spl1[0] } = $spl1[1]; ##map transcript to ensembl ID

    if ( (exists $geneToChrom{ $spl1[1] }) && ( $geneToChrom{ $spl1[1] } ne $spl1[3]) ) {
	print "Problem with ".$spl1[1]."\n";
    } else { 
	$geneToChrom{ $spl1[1] } = $spl1[3];
    }

    $geneToStrand{ $spl1[1] } = $spl1[4];
    $transcriptToChrom{ $spl1[0] } = $spl1[3];

}
close (INP1);
#print $transcriptToChrom{ "XM_001472958" }, "\n"; exit;
##################################################

my %problem = ();

my $logfile = "logfile_gtf_perl_".$species.".out";

print "Reading the file $inputGTF\n";

open (LOG, " > $logfile") or die;
open (OUT, " > $output ") or die;
open (INP2, " < $inputGTF ") or die;

print "Log in $logfile\n";

while (<INP2>) {
    #print $_;
    my @spl = split(' ', $_);

    my $transcript = $spl[ $transcriptNB ];
    $transcript =~ s/\"//g;
    $transcript =~ s/\.[0-9];//g;


    if ($transcript eq "NM_001042376") {
	print LOG "$transcript in the INS region, needs to be removed\n";
	next;
    }
    
    #print "Looking for $transcript\n"; exit;
    if (exists $hash{ $transcript }) {
	#print $transcript;
	if ($transcriptToChrom{ $transcript } ne $spl[0]) {
	    if (!exists $problem {$transcript }) {print LOG "Problem with $transcript which biomart puts on chromosome ".$transcriptToChrom{ $transcript }." but the gff sees it on chromosome ".$spl[0]."\n";}
	    $problem {$transcript } = "yes";
	} else {

	    my $oldgene = $spl[ 9 ];
	    $oldgene =~ s/\"|;//g;
	    my $gene = $hash{ $transcript };
	    
	    if ( $geneToStrand{ $gene } ne $spl[6] ) {
		if (!exists $problem {$transcript }) {print LOG "Strand problem for $gene, annotated by Biomart on ".$geneToStrand{ $gene }." but the gff sees it on ".$spl[6]."\n";}
		$problem {$transcript } = "yes";
	    } else {
		$_ =~ s/\"$oldgene\"/\"$gene\"/g;
		print OUT $_;
	    }
	}
    }

}


close(INP2);
close (OUT);
close (LOG);

print "Output in $output\n";
