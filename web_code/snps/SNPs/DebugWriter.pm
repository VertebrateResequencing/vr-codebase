#
# Author:    Petr Danecek (pd3@sanger.ac.uk)    Team 145
#
#--------------- DebugWriter -------------------------------------------
#
# The same as Writer, except no header nor footer is produced. This is used
#   to check the HTML code validity. (The Sanger code does not pass the validation.)
#

package SNPs::DebugWriter;

use strict;
use warnings;
use base qw(SNPs::Writer);

sub new
{
    my ($class, $args) = @_;
    my $self = $class->SUPER::new($args);
    if ( !exists($$self{'fname'}) ) { $$self{'fname'}=''; }
    return $self;
}

sub header
{
    my ($self,@args) = @_;

    my $css_styles='';
    for my $cssfile (@{$$self{'css'}})
    {
        $css_styles .= qq[<link rel="stylesheet" type="text/css" href="$cssfile" />\n];
    }

    my $js_scripts='';
    for my $jsfile (@{$$self{'jsfile'}})
    {
        $js_scripts .= qq[<script type="text/javascript" src="$jsfile" ></script>\n];
    }

    return qq[Content-type: text/html; charset=utf-8

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<title>Hello</title>
$css_styles
$js_scripts
</head>
<body>
    ];
}

sub print_footer
{
    my ($self,@args) = @_;
    $self->out("</body></html>");
}

1;


