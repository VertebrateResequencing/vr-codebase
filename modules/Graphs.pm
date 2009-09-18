package Graphs;

use strict;
use warnings;
use Utils;

our $R_CMD = '/software/R-2.9.0/bin/R';

=pod

=head1 METHODS

=head2 plot_stats

        Description: Create graphs for the given data sets.
        Arg [1]    : The statistics collected by e.g. SamTools::collect_detailed_bam_stats. 
                     Multiple graphs can be produced in a single run. The following data 
                     structure is expected:
                            outfile         .. file types supported by R: png, pdf, or jpg
                            title           .. 
                            desc_xvals      .. the x-axis label 
                            desc_yvals      .. the y-axis label
                            barplot         .. 1 for barplot
                            normalize       .. if set, yvals will be scaled so that the max value is 1
                            r_cmd           .. extra R commands to be executed
                            r_plot          .. extra R statements to plot(), such as e.g. "xlim=c(0,10)"
                            data => \[      .. Multiple lines can be plotted in one graph
                                { xvals=>xvals1, yvals=>yvals1 },
                                { xvals=>xvals2, yvals=>yvals2 },
                                ...
                            ],
        Returntype : None

=cut

sub plot_stats
{
    my ($stats) = @_;

    if ( !exists($$stats{'data'}) ) { Utils::error("FIXME\n") }
    if ( !exists($$stats{'outfile'}) ) { Utils::error("FIXME\n") }

    if ( !($$stats{'outfile'} =~ /([^.]+)$/ ) ) # Match the suffix
    { 
        Utils::error("Could not determine the filetype of \"$$stats{'outfile'}\".\n") 
    }
    my $file_type = $1;

    my $title   = exists($$stats{'title'}) ? $$stats{'title'} : '';
    my $xlabel  = exists($$stats{'desc_xvals'}) ? $$stats{'desc_xvals'} : '';
    my $ylabel  = exists($$stats{'desc_yvals'}) ? $$stats{'desc_yvals'} : '';
    my $barplot = exists($$stats{'barplot'}) ? 1 : 0;

    if ( $barplot && @{$$stats{'data'}} > 1 ) { Utils::error("TODO: multiple barplots in one graph.\n"); }

    open(my $fh, '>', "$$stats{'outfile'}.R") or Utils::error("$$stats{'outfile'}.R: $!");

    print $fh "par(bg='cornsilk')\n";
    print $fh "$file_type(file='$$stats{'outfile'}')\n";

    my $set  = 0;
    for my $vals (@{$$stats{'data'}})
    {
        if ( !$vals ) { Utils::error("Given ampty data set for $$stats{'outfile'}\n") }
        if ( !scalar @{$$vals{'xvals'}} ) { Utils::error("The data set is empty for $$stats{'outfile'}\n") }

        my $xrange = "x$set";
        my $yrange = "y$set";
        if ( $set > 0 )
        {
            $xrange = "xrange,x$set";
            $yrange = "yrange,y$set";
        }

        my ($x,$y);
        if ( $barplot )
        {
            $x = "'" . join("','", @{$$vals{'xvals'}}) . "'";
        }
        else
        {
            $x = join(',', @{$$vals{'xvals'}});
        }

        if ( $$stats{normalize} )
        {
            my ($extreme);
            $extreme = $$vals{yvals}->[0];
            for my $y (@{$$vals{'yvals'}})
            {
                if ( abs($y)>$extreme ) { $extreme=abs($y); }
            }

            my @scaled = ();
            for my $y (@{$$vals{'yvals'}})
            {
                push @scaled, 1.0*$y/$extreme;
            }
            $y = join(',', @scaled);
        }
        else { $y = join(',', @{$$vals{'yvals'}}); }

        print $fh qq[
x$set <- c($x)
y$set <- c($y)

xrange <- range($xrange)
yrange <- range($yrange)

];
        $set++;
    }

    print $fh "par(cex=1.25)\n";
    if ( $barplot )
    {
        print $fh "par(las=2)\n";
        print $fh "barplot(y0,,,x0, xlab='$xlabel',ylab='$ylabel')\n";
    }
    else
    {
        my $r_plot = exists($$stats{'r_plot'}) ? ", $$stats{'r_plot'}" : '';
        print $fh "plot(xrange,yrange,type='n', xlab='$xlabel',ylab='$ylabel' $r_plot)\n";
        $set = 0;
        for my $vals (@{$$stats{'data'}})
        {
            my $lines = exists($$vals{lines}) ? $$vals{lines} : '';
            my $type  = exists($$vals{type}) ? "type='$$vals{type}'" : "type='l'";
            print $fh "lines(x$set,y$set,$type,col=",$set+1,",pch=",$set+1,"$lines)\n";
            $set++;
        }
    }
    print $fh "title(main='$title', font.main=3)\n";
    print $fh $$stats{'r_cmd'} unless !exists($$stats{'r_cmd'});

    close $fh;

    Utils::CMD("cat $$stats{'outfile'}.R | $R_CMD --slave --vanilla");
    return;
}


sub create_gc_depth_graph
{
    my ($bindepth_file,$gcdepth_R,$png_file) = @_;

    # Parse the .bindepth file and create stats so that the R script can generate the graph
    #   - basically counts the lines and determines the bin size.
    #
    open(my $fh,'<',$bindepth_file) or Utils::error("$bindepth_file: $!");
    my $nlines   = 0;
    my ($prev_pos,$pos,$bin_size,$chrm,$chrm_prev);
    while (my $line=<$fh>)
    {
        $nlines++;

        # 1       10000   0       NA      0.0000  1.0000
        if ( !($line=~/^(\S+)\s+(\d+)\s+/) ) { Utils::error("Expected different output: $line"); }

        $chrm = $1;
        $pos  = $2;
        if ( defined $prev_pos && $chrm eq '1' )
        {
            if ( !defined $bin_size ) { $bin_size = $pos - $prev_pos; }
            if ( $bin_size != $pos-$prev_pos )
            {
                Utils::error("The bin_size of diffent size on line $nlines: $bin_size vs ".($pos-$prev_pos)."\n");
            }
        }

        $prev_pos  = $pos;
        $chrm_prev = $chrm;
    }
    close $fh;

    # Finally, create the R script and run the command.
    open($fh, '>', "$png_file.R") or Utils::error("$png_file.R: $!");
    print $fh qq[
source('$gcdepth_R')
depdat = read.depdat('$bindepth_file', Ndat = $nlines, bin = $bin_size)
gcdepth(depdat, sname = '', depmax = NULL, hc = TRUE, nbins = 30, binned = TRUE, outfile = '$png_file')
];

    close $fh;
    Utils::CMD("cat $png_file.R | $R_CMD --slave --vanilla");

    return;
}


1;


