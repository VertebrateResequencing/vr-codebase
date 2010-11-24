#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most tests => 18;
    
    use_ok('VertRes::Parser::dict');
}

my $pd = VertRes::Parser::dict->new();
isa_ok $pd, 'VertRes::Parser::ParserI';
isa_ok $pd, 'VertRes::IO';
isa_ok $pd, 'VertRes::Base';

ok my $rh = $pd->result_holder(), 'result_holder returned something';
is ref($rh), 'HASH', 'result_holder returns a hash ref';
is keys %{$rh}, 0, 'the result_holder starts off empty';

ok ! $pd->next_result, 'next_result returns false when we have no file set';

my $d_file = File::Spec->catfile('t', 'data', 'human_b36_male.dict');
ok -e $d_file, 'file we will test with exists';
ok $pd->file($d_file), 'file set into parser';

ok $pd->next_result, 'next_result now works';
is_deeply $rh, {SN => 1, LN => 247249719, UR => 'file:/lustre/scratch103/sanger/team145/g1k/ref/human_b36_male.fa', M5 => '9ebc6df9496613f373e73396d5b3b6b6'}, 'result_holder contains correct info for first line';
ok $pd->next_result, 'next_result worked again';
is $rh->{M5}, 'b12c7373e3882120332983be99aeb18d', 'result_holder contains correct info for second line';

# check the last line as well
while ($pd->next_result) {
    next;
}
is_deeply $rh, {SN => 'NC_007605', LN => 171823, UR => 'file:/lustre/scratch103/sanger/team145/g1k/ref/human_b36_male.fa', M5 => '6743bd63b3ff2b5b8985d8933c53290a'}, 'result_holder contains correct qname for last line';

# seq_lengths
is_deeply {$pd->seq_lengths}, { qw(1	247249719
2	242951149
3	199501827
4	191273063
5	180857866
6	170899992
7	158821424
8	146274826
9	140273252
10	135374737
11	134452384
12	132349534
13	114142980
14	106368585
15	100338915
16	88827254
17	78774742
18	76117153
19	63811651
20	62435964
21	46944323
22	49691432
X	154913754
Y	57772954
MT	16571
NT_113887	3994
NT_113947	4262
NT_113903	12854
NT_113908	13036
NT_113940	19187
NT_113917	19840
NT_113963	24360
NT_113876	25994
NT_113950	28709
NT_113946	31181
NT_113920	35155
NT_113911	36148
NT_113907	37175
NT_113937	37443
NT_113941	37498
NT_113909	38914
NT_113921	39615
NT_113919	40524
NT_113960	40752
NT_113945	41001
NT_113879	42503
NT_113938	44580
NT_113928	44888
NT_113906	46082
NT_113904	50950
NT_113873	51825
NT_113966	68003
NT_113943	81310
NT_113914	90085
NT_113948	92689
NT_113886	96249
NT_113932	104388
NT_113929	105485
NT_113878	106433
NT_113927	111864
NT_113900	112804
NT_113918	113275
NT_113875	114056
NT_113942	117663
NT_113926	119514
NT_113934	120350
NT_113954	129889
NT_113953	131056
NT_113874	136815
NT_113883	137703
NT_113924	139260
NT_113933	142595
NT_113884	143068
NT_113890	143687
NT_113870	145186
NT_113881	146010
NT_113939	147354
NT_113956	150002
NT_113951	152296
NT_113902	153959
NT_113913	154740
NT_113958	158069
NT_113949	159169
NT_113889	161147
NT_113936	163628
NT_113957	166452
NT_113961	166566
NT_113925	168820
NT_113882	172475
NT_113916	173443
NT_113930	174588
NT_113955	178865
NT_113944	182567
NT_113901	182896
NT_113905	183161
NT_113872	183763
NT_113952	184355
NT_113912	185143
NT_113935	185449
NT_113880	185571
NT_113931	186078
NT_113923	186858
NT_113915	187035
NT_113885	189789
NT_113888	191469
NT_113871	197748
NT_113964	204131
NT_113877	208942
NT_113910	211638
NT_113962	217385
NT_113899	520332
NT_113965	1005289
NT_113898	1305230
NC_007605	171823) }, 'seq_lengths worked';

# try out a new dict file with more tags
$d_file = File::Spec->catfile('t', 'data', 'ncbi37.dict');
ok -e $d_file, 'file we will test with exists';
$pd->file($d_file);
$pd->next_result;
is_deeply $rh, {SN => 1, LN => 249250621, AS => 'NCBI37', UR => 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz', M5 => '1b22b98cdeb4a9304cb5d48026a85128', SP => 'Human'}, 'result_holder contains correct info for first line of an ncbi37 dict';
