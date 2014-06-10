#!/usr/local/bin/perl
# Generate dbf of external gene names for use by vcf2consequences_vep --n option

BEGIN{
    my $ROOT = '/software/vertres/installs';
    my $VERSION = '75';
    unshift(@INC, "$ROOT/bioperl-1.2.3/lib/site_perl/5.8.8");
    unshift(@INC, "$ROOT/ensembl/$VERSION/ensembl/modules");
    unshift(@INC, "$ROOT/ensembl/$VERSION/ensembl-variation/modules");
    unshift(@INC, "$ROOT/ensembl/$VERSION/ensembl-compara/modules");
    unshift(@INC, "$ROOT/ensembl/$VERSION/ensembl-functgenomics/modules");
}

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use DB_File;

die "Usage $0 <species> <db_filename>" unless @ARGV == 2;
my $species = $ARGV[0];
my $dbf = $ARGV[1];

my $reg = "Bio::EnsEMBL::Registry";
$reg->load_registry_from_db( -host => "ensembldb.ensembl.org", -user => "anonymous",);
my $ga = $reg->get_adaptor("$species", "core", "Gene");
my $sa = $reg->get_adaptor("$species", "core", "Slice");

my %TBL;
$DB_HASH->{'cachesize'}=2000000000;
$DB_HASH->{'nelem'}=100000; # estimate no of elements
tie %TBL, "DB_File", "$dbf", O_RDWR|O_CREAT, 0666, $DB_HASH or die "Cannot open file : $!\n";

my $c=0;
foreach my $chr_slice (@{$sa->fetch_all('toplevel')}){
   foreach my $gene (@{$ga->fetch_all_by_Slice($chr_slice)}){
		if (defined($gene->external_name())) {
            $TBL{$gene->stable_id} = $gene->external_name();
            $c++;
		}
   }
}
print STDERR  "Wrote $c external gene names for $species to $dbf\n";
