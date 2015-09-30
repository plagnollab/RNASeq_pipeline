#!/share/apps/perl-5.14.2/bin/perl

use lib '/cluster/project8/vyp/vincent/libraries/perl/ensembl_81/bioperl-live';
use lib '/cluster/project8/vyp/vincent/libraries/perl/ensembl_81/ensembl/modules';
use lib '/cluster/project8/vyp/vincent/libraries/perl/ensembl_81/ensembl-variation/modules';


use DBI;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;

use strict;

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
				 );

sub feature2string
{
    my $feature = shift;

    my $stable_id  = $feature->stable_id();
    my $seq_region = $feature->slice->seq_region_name();
    my $start      = $feature->start();
    my $end        = $feature->end();
    my $strand     = $feature->strand();

    return sprintf( "%s: %s:%d-%d (%+d)",
        $stable_id, $seq_region, $start, $end, $strand );
}	


my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
my $tr_adaptor    = $registry->get_adaptor( 'Human', 'Core', 'Transcript' );


my @list = ("X", "Y", (1..22));


foreach my $chr (@list) {
# Get a slice of chromosome 20, 10MB-11MB
    my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr, 0, 11e8 );
    


    my $genes = $slice->get_all_Genes();
    while ( my $gene = shift @{$genes} ) {
	my $transcripts = $gene->get_all_Transcripts();
	while ( my $transcript = shift @{$transcripts} ) {
	    my $strand = $transcript->strand();
	    my $cleanStrand = "NA";
	    if ($strand == -1){ $cleanStrand = "-";}
	    if ($strand == 1) {$cleanStrand = "+";}

	    print $chr."\t".$gene->stable_id()."\t".$transcript->stable_id()."\t".$cleanStrand."\n";
	}
    }
}
