#!/usr/bin/env perl
use strict;
use warnings;

# supply lane name and qc_status string
# we then set that status on the node for that lane in VRPipe's graph db
# we do this via VRPipe's webserver since loading VRPipe code is slow

my ($lane, $status) = @ARGV;
chomp($status);

system("wget -q --no-check-certificate https://vr-2-2-02.internal.sanger.ac.uk:9091/qcgrind_lane_status_update/$lane/$status -O /dev/null");

exit;
