#!/usr/local/bin/perl -T
#
# ~pd3/git/wscripts/bstats -c /nfs/WWWdev/INTWEB_docs/lib/teams/team145/vr-codebase/web_code/bstats/dat/config -d /nfs/WWWdev/INTWEB_docs/lib/teams/team145/vr-codebase/web_code/bstats/dat/uk10k
# ~pd3/git/wscripts/bstats -c /nfs/WWWdev/INTWEB_docs/lib/teams/team145/vr-codebase/web_code/bstats/dat/config -d /nfs/WWWdev/INTWEB_docs/lib/teams/team145/vr-codebase/web_code/bstats/dat/farm2
#
#
use strict;
use warnings;
use URI;
use CGI;

use SangerPaths qw(core team145);
use SangerWeb;
use VertRes::QCGrind::ViewUtil;

$ENV{PATH}="";

my $opts = 
{
    # Real path is /nfs/WWWdev/INTWEB_docs/lib/teams/team145/vr-codebase/web_code/bstats
    config  => 'dat/config',
    wwwdir  => '/Teams/Team145/bstats',
};

my $utl = VertRes::QCGrind::ViewUtil->new();
my $title = 'bstats';
my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'style'   => $utl->{CSS},
});
   
$$opts{cgi} = $sw->cgi();

print $sw->header();
print qq[<div class="centerFieldset"'>
       ];

main($opts);

print qq[
    </div>];
   
print $sw->footer();

exit;

#---------------------------------------

sub main
{
    my ($opts) = @_;
    
    my $list_users = read_config($$opts{config});
    my $list_farms = { 'uk10k'=>1, 'farm2'=>1 };
    my $list_views = { day=>1, week=>1, month=>1, all=>1 };
    my $q = $$opts{cgi};
    my $user = $q->param('u');
    my $farm = $q->param('f');
    my $view = $q->param('v');
    if ( !$user or !exists($$list_users{$user}) ) { $user='all' } 
    if ( !$farm or !exists($$list_farms{$farm}) ) { $farm='uk10k' } 
    if ( !$view or !exists($$list_views{$view}) ) { $view='day' }
    $user = ($user=~/([\w-]+)/) ? $1 : '';
    $farm = ($farm=~/([\w-]+)/) ? $1 : '';
    $view = ($view=~/([\w-]+)/) ? $1 : '';
    my $header = qq[
                <form action="" method="post" style="padding-bottom:0.5em;">
                User: <select name="u" style="border:0px solid;margin:0px;padding:0px;" id="users">
        ];
    for my $usr (sort keys %$list_users)
    {
        my $selected = $usr eq $user ? ' selected' : '';
        $header .= qq[<option$selected>$usr</option>\n];
    }
    $header .= qq[
                </select>
                <span style="margin-left:1em;">Farm:</span> <select name="f" style="border:0px solid;margin:0px;padding:0px;" id="farms">
            ];
    for my $frm (sort keys %$list_farms)
    {
        my $selected = $frm eq $farm ? ' selected' : '';
        $header .= qq[<option$selected>$frm</option>\n];
    }
    $header .= qq[
                </select>
                </form>
        ];
    my @titles;
    my @imgs;
    my $uri;
    for my $file qw(day week month all)
    {
        if ( !-e "dat/$farm/$file.$user.png" ) { next; }
        $uri = URI->new("data:");
        $uri->media_type("image/png");
        $uri->data(scalar(`/bin/cat dat/$farm/$file.$user.png`));
        push @titles, $file;
        push @imgs, qq[
                    <img 
                        src="$uri" 
                        style="border:1px solid;cursor:pointer;width:100px;" 
                        onmouseover="
                            current_view = '$file';
                            \$('#main_img').attr('src','$uri');"
                        >
                 ];
    }

    my $ncols = scalar @titles;
    print qq[
        <table>
        <tr><td>
                $header 
            </td></tr>
        <tr><td>
                <a href="http://www.wtgc.org/cgi-bin/IT/ISG/diskusage?server=new_lustre">Lustre limits</a>
                <br>
                <a href="http://farm2-srv1.internal.sanger.ac.uk/rrd/cgi/lustre.cgi?_submitted=68&timespan=1hour&filesystem=lus02&cfactor=AVERAGE">Lustre load</a>
            </td></tr>
        <tr><td>
            <table style="border-collapse:collapse;margin-left:auto;margin-right:auto;padding:1em;background-color:white;">
            <tr class="titles"><td>] .join('</td><td>',@titles). q[</td></tr>
            <tr class="imgs"><td>] .join('</td><td>',@imgs). qq[</td></tr>
            </table>
            </td></tr>
        <tr><td>
                <img src="$uri" id="main_img">
            </td></tr>
        </table>

        ] . q[

        <script src="http://code.jquery.com/jquery-latest.min.js" type="text/javascript"></script>
        <script type="text/javascript">
            function getLink()
            {
                return '?u=' + $("#users option:selected").text() + '&f=' + $("#farms option:selected").text() + '&v=' + current_view;
            }
            $('#users').change(function() {
                window.location.href = getLink();
            })
            $('#farms').change(function() {
                window.location.href = getLink();
            })
            current_view = '] . $view . qq[';
        </script>
    ];
}

sub read_config
{
    my ($file,$regex) = @_;
    my %out;
    open(my $fh,'<',$file) or die "$file: $!";
    while (my $line=<$fh>)
    {
        chomp($line);
        $line = ($line=~/([\w-]+)/) ? $1 : '';
        $out{$line} = 1;
    }
    close($fh);
    return \%out;
}

