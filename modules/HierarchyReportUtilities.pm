package HierarchyReportUtilities;
use strict;
use Carp;
use File::Basename;
use File::Spec; 
use File::Copy;
use Cwd;

=pod

=head1 NAME

=head1 SYNOPSIS

a function to produce a summary of one column vs. another

=head1 ARGUMENTS

=head1 DESCRPTION

=head1 AUTHOR

Thomas Keane, I<tk2@sanger.ac.uk>

=cut
sub summariseValueByColumn
{
	croak "Usage: summariseValueByColumn summary_csv selectColumn summaryColumn valueColumn operation output R_script" unless @_ == 6;
	my $csv = shift;
	my $selectColumn = shift;
	my $summaryCol = shift;
	my $valueCol = shift;
	my $operation = shift;
	my $output = shift;
	my $rscript = shift;
	
	croak "Cant find csv file: $csv\n" unless -f $csv;
	croak "Invalid operation: $operation\n" unless $operation eq 'plus';
	
	#build a hash of arrays
	my %summary;
	open( CSV, $csv ) or die $!;
	<CSV>; #header
	while( <CSV> )
	{
		chomp;
		my @s = split( /,/, $_ );
		
		if( $s[ $selectColumn ] 
		if( defined $summary{ $s[ $summaryCol ] } )
		{
			if( $operation eq 'plus' )
			{
				$summary{ $s[ $summaryCol ] } += $s[ $valueCol ];
			}
		}
		else
		{
			$summary{ $s[ $summaryCol ] } = $s[ $valueCol ];
		}
	}
	close( CSV );
	
	open( OUT, ">$output" ) or die "Cannot craete output: $!\n";
	foreach( keys( %summary ) )
	{
		print "$_ $summary{ $_ }\n";
	}
	close( OUT );
}

1;
