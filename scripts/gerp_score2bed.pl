#!/usr/local/bin/perl
# generate GERP scores, write to bed file
# accepts species amd chr as parameters
BEGIN{
    my $ROOT = '/software/vertres/lib/all';
    my $VERSION = '70';
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
use Data::Dumper;

die "Usage $0 <species> <chromosome>" unless @ARGV == 2;
my $species = $ARGV[0];
my $chr = $ARGV[1];
print STDERR `date`, "Generating bed for $species, chr $chr\n";

my $reg = "Bio::EnsEMBL::Registry";
$reg->load_registry_from_db( -host => "ensembldb.ensembl.org", -user => "anonymous",);

my $mlss_adaptor = $reg->get_adaptor("Multi", "compara", "MethodLinkSpeciesSet");
my $mlss = $mlss_adaptor->fetch_by_method_link_type_species_set_name("GERP_CONSERVATION_SCORE", "mammals");
my $cs_adaptor = $reg->get_adaptor("Multi", 'compara', 'ConservationScore');
my $va_adaptor = $reg->get_adaptor("$species", 'variation', 'variation');
my $vf_adaptor = $reg->get_adaptor("$species", 'variation', 'variationfeature');
my $slice_adaptor = $reg->get_adaptor("$species", 'core', 'slice');

my $gerp_slice = $slice_adaptor->fetch_by_region('chromosome', $chr);

my $slice_len =  $gerp_slice->end;
print STDERR "$chr:($slice_len)\n";

my $slice_sz = 100000;
my $i;
for ($i=1;$i<$slice_len;$i+=$slice_sz) {
    my $end = $i+$slice_sz-1;
    $end = $slice_len if $end > $slice_len;
    $gerp_slice = $slice_adaptor->fetch_by_region('chromosome', $chr,$i,$end);
    my $offset = $gerp_slice->start-1;

    #To get one score per base in the slice, must set display_size to the size of the slice.
    my $display_size = $end - $i;
    my $scores = $cs_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $gerp_slice, $display_size);

    print STDERR `date`, "$chr:Number of scores ($i - $end, $display_size) = " . @$scores . "\n" if @$scores > 0;

    foreach my $score (@$scores) {
        if (defined $score->diff_score) {
            printf("$chr\t%d\t%d\t%.4f\n",  $score->position + $offset, $score->position + $offset, $score->diff_score);
        }
    }
}
