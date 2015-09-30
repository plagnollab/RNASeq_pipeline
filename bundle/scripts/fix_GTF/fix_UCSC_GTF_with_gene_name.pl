#!/share/apps/perl-5.14.2/bin/perl

#use strict;   

my $species = $ARGV[0];

my $biomart;
my $output;
my $inputGTF;
my $transcriptNB = 11;

if ($species eq "human_hg38") {
    $biomart = "human_hg38/transcript_gene_info.tab";
    $output = "/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GTF_files/human_hg38_fixed.gtf";
    $inputGTF = "/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GTF_files/human_hg38.gtf";
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
    $hash{ $spl1[2] } = $spl1[1]; ##map transcript to ensembl ID
    $geneToChrom{ $spl1[1] } = $spl1[0]; ##map gene to chromosome
    $geneToStrand{ $spl1[1] } = $spl1[3]; ##map gene to strand
}
close (INP1);

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
    $transcript =~ s/\"|;//g;
    $transcript =~ s/\.[0-9];//g;
    ##print $transcript."\n";

    if (exists $hash{ $transcript }) {
	my $oldgene = $spl[ 9 ];
	my $gene = $hash{ $transcript };
	my $chromosome = $geneToChrom{ $gene };
	my $strand = $geneToStrand{ $gene };
	
	my $gtfChrom = $spl[ 0 ];
	$gtfChrom =~ s/chr//g;
	
	my $gtfStrand = $spl[6];

	if ($chromosome ne $gtfChrom) {
	    print LOG "Problem with gene $gene that is not on chromosome $chromosome but instead on $gtfChrom.\n";
	} else {
	    if ($strand ne $gtfStrand) {
		print LOG "Problem with gene $gene that is not on strand $strand but instead on $gtfStrand\n";
	    } else {
		$_ =~ s/gene_id \"$transcript\"/gene_id \"$gene\"/g;
		print OUT $_;
	    }
	}
    }

}


close(INP2);
close (OUT);
close (LOG);

print "Output in $output\n";
