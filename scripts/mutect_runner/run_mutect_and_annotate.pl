#!/usr/bin/env perl

# Author: Kim Wong kw10@sanger.ac.uk

use strict;
use warnings;
use Carp;
use Utils;
use FindBin qw($Bin);

my $bin = $Bin;

# input required: config file
# ./run_mutect_and_annotate.pl +sampleconf
# to get an example
# To run:
# run_mutect_and_annotate.pl +config mutect.conf -o $PWD +loop 150 +mail kw10@sanger.ac.uk +nocache &>log
# Requires scripts from the mutect_runner package:
# reformatVCF.pl
# check_samples.pl

# Create new runner object and run it
my $runner = mutectRunner->new();
$runner->run();

##exit;

#-----------------------

package mutectRunner;
use base qw(Runner);
use strict;
use warnings;
use Getopt::Std;


sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    $$self{do_clean} = 0;
    $$self{limits} = { memory=>5_000, runtime=>24*60 };
    $$self{mutect} = '/software/vertres/bin-external/mutect-1.1.7.jar';
	$$self{java} = '/software/vertres/bin-external/java7';
	#$$self{vep}    = '/software/vertres/bin-external/variant_effect_predictor.pl';
	##$$self{vep}    = '/nfs/users/nfs_k/kw10/src/variant_effect_predictor_ens75.pl';
	$$self{vep}    = '/software/vertres/bin-external/variant_effect_predictor_ens78.pl';
	$$self{bcftools} = '/software/vertres/bin-external/bcftools';
	$$self{reformat_vcf} = "$bin/reformatVCF.pl";
	$$self{mutect_summary} = "$bin/check_samples.pl";
    $$self{_sampleconf} = q[


			# This version of MuTect runner is compatible with MuTect 1.1.7

			# REQUIRED: paths to required software:
			# path to required runner scripts
			reformat_vcf =>  '] .$$self{reformat_vcf} .  q[',
			mutect_summary => '] . $$self{mutect_summary} . q[',
			# Note java 7 is required for muTect-1.1.7
			java  => '] .$$self{java}.  q[',
			mutect  => '] .$$self{mutect}.  q[',
			bcftools => '] .$$self{bcftools}. q[',
			# Ensure the VEP version matches the VEP cache and Ensembl API
			vep => '] .$$self{vep}. q[',
			# OPTIONAL: path to bedtools if using 'targets' key. Default is 'bedtools'
			bedtools => 'bedtools',
			# REQUIRED: tab-delimited file of 
			# tumSampleName,NormalSampleName,tumourBamPath,normalBamPath
			bams     => 'bams.list',

			# REQUIRED: 'y' or 'n'; run mutect faster by running each chrom in parallel
			# Highly recommended for WGS data since MuTect can be slow
			bychrom => 'y',

			# OPTIONAL: mutect options other than cosmic file, snp file, wig file, chromosome, reference file, input/output
			mutect_opt => "--max_alt_allele_in_normal_fraction 0.1 --pir_median_threshold 5 --heavily_clipped_read_fraction 0.25",

			## Comment/alter or delete this section if not running human data ----##
			## VEP options, example for human:
			# REQUIRED: reference genome
			reference   => '/lustre/scratch105/projects/g1k/ref/main_project/human_g1k_v37.fasta',
			# OPTIONAL: dbSNP vcf and Cosmic VCF for MuTect
			# Note: using these may affect whether a call is filtered or passed; they are not used simply
			# for annotating SNV calls. See MuTect website for details. Update files as required.
			#snpfile    =>  '/lustre/scratch105/vrpipe/refs/human/ncbi37/dbsnp/dbsnp_141.b37.vcf.gz',
			#cosmic   => '/lustre/scratch105/vrpipe/refs/human/ncbi37/resources/Cosmic/20140805/CosmicCoding_and_NonCoding.vcf.gz',
			# REQUIRED: VEP parameters
			# Note that up to Ensembl 75 (Feb 2014) is GRCh37.p13; starting at Ensembl 76 is GRCh38
			# The Ensemble vep_cache version must match the api version
			ens_version => '75',
			species => 'human',
			assembly => 'GRCh37',
			ensembl_api => '/software/vertres/installs/ensembl/75/',
			vep_cache => '/lustre/scratch105/projects/g1k/ref/vep_cache',
			vep_fasta => undef,     # optional --fasta argument for variant_effect_predictor.pl
			vep_stats => 'n',
			add_hgvs => 'y',
			chroms => [ qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y) ],
			##--------------------------------------------------------------------##

			## Comment or delete this section if not running mouse data ----------##
			## VEP options, example for mouse:
			# REQUIRED: reference genome
			reference   => '/lustre/scratch105/vrpipe/refs/mouse/GRCm38/GRCm38_68.fa',
			# OPTIONAL: dbSNP vcf for MuTect
			# Note: using a SNP file may affect whether a call is filtered or passed; they are not used simply
			# for annotating SNV calls. See MuTect website for details. Update SNP file as required.
			#snpfile    =>  '/lustre/scratch105/vrpipe/refs/mouse/GRCm38/dbSNP/138/all.snps.dbSNP138.vcf.gz',
			# REQUIRED: VEP parameters
			ens_version => '78',
			species => 'mouse',
			ensembl_api => '/software/vertres/installs/ensembl/78/',
			vep_cache => '/lustre/scratch105/projects/g1k/ref/vep_cache',
			vep_stats => 'n',
			chroms => [ qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y) ],
			##--------------------------------------------------------------------##


			# Output filtering and formatting
			# REQUIRED: 'y' or 'n'
			# DP4: allele counts per strand for the tumour
			add_DP4T  => 'y',
			
			# REQUIRED: 'y' or 'n'
			# Remove MuTect REJECT calls
			remove_fail => 'y',
			
			# OPTIONAL: BED file for target region, eg, for exomes
			# targets => '', 
			
			# REQUIRED: 'vcf' or 'table' or 'both'
			# Output a vcf/tab file for each tumour/normal MuTect run
			output => 'both',
			
			# REQUIRED" 'y' or 'n'
			# For mutlitple tum/norm MuTect runs, create an optional summary table of calls 
			# from all tumour/norm pairs. This is a table listing all sites of interest and 
			# annotations with one column per tumour and 1 or 0 indicating presence/absence
			summary_table => 'y',
			
			# REQUIRED: 'default', 'all' or a list of "|" separated SO consequences used in VEP
			# http://www.ensembl.org/info/genome/variation/predicted_data.html#consequences
			# 'default' is 'transcript_ablation|stop_g|stop_l|missense|splice_d|splice_a|initiator_c|start_lost|transcript_amp|frameshift_v|inframe_'
			conseq_list    => 'default',
			# OPTIONAL: filter by transcript biotype
			# Keep consequences with transcripts of these Ensembl biotypes; separate with "|"
			# biotype_filter => 'protein_coding',
	]."\n";


	

    $$self{usage} .=
        "Usage: run-mutect_and_vep\n" .
        "Options:\n" .
        "   -c, --clean             Clean all temporary files (and do nothing else)\n" .
        "   -o, --outdir <dir>      Output directory [required]\n" .
        "\n";

    return $self;
	

}
# from config file:

sub parse_args {
	my ($self) = @_;
	my @required = (
		'add_DP4T',
		'add_hgvs',
		'bams',
		'bcftools',
		'bychrom',
		'chroms',
		'conseq_list',
		'ensembl_api',
		'ens_version',
		'java',
		'mutect',
		'output',
		'reference',
		'remove_fail',
		'species',
		'summary_table',
		'vep',
		'vep_cache',
		'vep_stats',
	);
	my @choose = (
		'bychrom',
		'vep_stats',
		'add_DP4T',
		'remove_fail',
		'summary_table',
		'add_hgvs',
	);
	my @files = (
		'bams',
		'bcftools',
		'ensembl_api',
		'java',
		'mutect',
		'reference',
		'vep',
		'vep_cache',
	);
	my @optional = (
		'mutect_opt',
		'targets',
		'snpfile',
		'biotype_filter',
		'cosmic',
		'bedtools',
		'assembly'
	);

    while (defined(my $arg=shift(@ARGV))) {
        if ( $arg eq '-o' or $arg eq '--outdir' ) {
			$$self{outdir}=shift(@ARGV);
			next;
		}
		$self->throw();
    }
    if ( !exists($$self{outdir}) ) { $self->throw("Missing the -o option.\n"); }
	# check for require parameters in config file:
	foreach my $op (@required) {
		if (!exists($$self{$op})) {
			$self->throw("Required key '$op' in config file is missing");
		}
		if ($op eq 'output') {
			$$self{$op} = lc($$self{$op});
			if ($$self{$op} ne 'both' && $$self{$op} ne 'table' && $$self{$op} ne 'vcf') {
				$self->throw("Option 'output' must be 'both' or 'table' or 'vcf'");
			}
		}
	}
	foreach my $c (@choose) {
		$$self{$c} = lc($$self{$c});
		$$self{$c} =~ s/yes/y/;
		$$self{$c} =~ s/no/n/;
		if ($$self{$c} ne 'y' && $$self{$c} ne 'n') {
			$self->throw("Option $c requires y(es) or n(o) in config file");
		}
	} 
	foreach my $f (@files) {
		if (! -e $$self{$f}) {
			$self->throw("No such file $$self{$f}");
		}
	}
	if ($$self{targets} && !$$self{bedtools}) {
		$$self{bedtools} = 'bedtools';
	}
	foreach my $o (@optional) {
		if ($$self{$o} && $o ne 'biotype_filter' && $o ne 'bedtools' && $o ne 'mutect_opt' && $o ne 'assembly') {
			if ( ! -e $$self{$o} ) {
				$self->throw("No such file $$self{$o}");
			}
		}
		if ( exists $$self{$o} && defined $$self{$o} ) {
			print STDERR "Optional setting: $o $$self{$o}\n"
		}
		else {
			print STDERR "Optional setting: $o NOT SET\n"
		}
	}
	if ( ! -e "$$self{vep_cache}/$$self{species}" ) {
		$self->throw("No directory $$self{species} in vep cache dir $$self{vep_cache}");
	}

}


sub main {
    my ($self) = @_;
    $self->parse_args();
    $self->save_config();

    my $outdir = $$self{outdir};
    my @chroms = @{$$self{chroms}};

    if ( $self->is_finished("$outdir/all_done") ) { $self->all_done; }

	$self->set_limits(%{$$self{limits}});
	my @mutect_list = `cut -f 1,2 $$self{bams} | sed 's/\\t/-vs-/'`;
	chomp @mutect_list;

	# run mutect on all tumour/normal pairs
	open F, "<$$self{bams}" || die "Cannot open bam file $$self{bam}\n";
	my %samples;
	while (<F>) {	
		my ($tname,$nname,$tbam,$nbam) = split "\t", $_;
		if (!$tname || !$nname || !$tbam || !$nbam) {
			$$self->throw("Malformed bam file list. Format is tumourname[tab]normalname[tab]tumourbam[tab]normalbam\n");
		}
		$samples{"$tname-vs-$nname"}{tbam}=$tbam;
		$samples{"$tname-vs-$nname"}{nbam}=$nbam;
	}

	# run mutect by chromosome and then merge files
	if ($$self{bychrom} eq 'y') {
		foreach my $sample (sort keys %samples) {
			`mkdir -p $outdir/$sample`;
			my $prefix = "$outdir/$sample/$sample";
			for my $chr (@{$$self{chroms}}) {
				$self->spawn('mutect_call',"$outdir/$sample/mutect.$chr.done",$prefix,$samples{$sample}{tbam},$samples{$sample}{nbam},$chr);
			}
		}
	}
	# run on whole bam:
	else {
		foreach my $sample (sort keys %samples) {
			`mkdir -p $outdir/$sample`;
			my $prefix = "$outdir/$sample/$sample";
			$self->spawn('mutect_call',"$outdir/$sample/mutect.done",$prefix,$samples{$sample}{tbam},$samples{$sample}{nbam});
		}
	}
	$self->wait;
	# merge the mutect output (stats and vcfs) if run by chrom
	if ($$self{bychrom} eq 'y') {
		foreach my $sample (sort keys %samples) {
			$self->spawn('mutect_merge_vcf',"$outdir/$sample/$sample.vcf.gz","$outdir/$sample",$sample);
			$self->spawn('mutect_merge_stats',"$outdir/$sample/$sample.mutect_stats.txt","$outdir/$sample",$sample);
		}
	}
	$self->wait;

	##exit;


	# optional: parse out DP4 info from mutect text output and add to vcf; remove reject calls, optional
    $self->set_limits( queue=>'small', memory=>1_000, runtime=>undef );
	if ($$self{add_DP4T} eq 'y') {
		foreach my $sample (sort keys %samples) {
			my $out = "$outdir/$sample/$sample";
			my ($statfile,$vcfin,$vcfout) = $$self{bychrom} eq 'no' ? ("$out.mutect_stats.txt","$out.vcf","$out.filt.vcf.gz") : ("$out.mutect_stats.txt","$out.vcf.gz","$out.filt.vcf.gz");
		##	print STDOUT qq(FOUND: 'add_dp4',"$vcfout.gz",$$self{remove_fail},$statfile,$vcfin,$out)."\n\n";
			$self->spawn('add_dp4',"$vcfout",$$self{remove_fail},$statfile,$vcfin,$outdir,$sample);
		}
	}
	# get PASS calls only from vcf
	elsif ($$self{remove_fail} eq 'y') {
		foreach my $sample (sort keys %samples) {
			my $out = "$outdir/$sample/$sample";
			$self->spawn('filter_reject',"$out.filt.vcf.gz","$out.vcf");
		}
	}
	else {
		for (@mutect_list) {
			my $cmd = "cd $$self{outdir}/$_ && ln -sv $_.vcf $_.filt.vcf && bgzip -c $_.filt.vcf > $_.filt.vcf.gz && tabix -p vcf $_.filt.vcf.gz && cd $$self{outdir}";
			#$self->cmd($cmd);
			system($cmd) == 0 || die "$cmd failed : $!";
		}
	}
	$self->wait;

	# optional : filter out target region SNVs
	if ($$self{targets}) {
		foreach my $sample (sort keys %samples) {
			my $out = "$outdir/$sample/$sample";
			$self->spawn('filter_targets',"$out.filt.targets.vcf.gz","$out.filt.vcf.gz");
		}
	}
	$self->wait if $$self{targets};

	foreach my $dir (@mutect_list) {
		if ($$self{targets}) {
			$self->spawn('symlinktarg',"$$self{outdir}/$dir/all.final.vcf.gz",$dir);
			##my $cmd = "cd $$self{outdir}/$dir && ln -svf $dir.filt.targets.vcf.gz all.final.vcf.gz && ln -svf $dir.filt.targets.vcf.gz.tbi all.final.vcf.gz.tbi && cd $$self{outdir}";
			##system($cmd) == 0 || die "$cmd failed : $!";
		}
		else {
			###print STDOUT "FOUND: cd to $$self{outdir}/$dir\n";
			$self->spawn('symlink',"$$self{outdir}/$dir/all.final.vcf.gz",$dir);
			##my $cmd = "cd $$self{outdir}/$dir && ln -svf $dir.filt.vcf.gz all.final.vcf.gz && ln -svf $dir.filt.vcf.gz.tbi all.final.vcf.gz.tbi && cd $$self{outdir}";
			##system($cmd) == 0 || die "$cmd failed : $!";
		}	
		
	}
	$self->wait;
	# for each tum/norm pair, run VEP by chromosome
	foreach my $sample (sort keys %samples) {
		my $out = "$outdir/$sample";
        for my $chr (@{$$self{chroms}}) {
			$self->spawn('split_vcf',"$out/all.$chr.vcf.gz","$out/all.final.vcf.gz",$chr);
		}
	}
	$self->wait;
    $self->set_limits( queue=>'normal', memory=>3_000, runtime=>undef );
	foreach my $sample (sort keys %samples) {
		my $out = "$outdir/$sample";
        for my $chr (@{$$self{chroms}}) {
			$self->spawn('run_vep',"$out/vep.$chr.done","$out/all.$chr.cons.vcf","$out/all.$chr.vcf.gz",);
		}
	}
	$self->wait;
	# check for missing chrom files (eg: if no chrY snps, vep does not output a file!
	# make bgzip and tabix
    $self->set_limits( queue=>'small', memory=>1_000, runtime=>undef );
	foreach my $sample (sort keys %samples) {
		my $out = "$outdir/$sample";
        for my $chr (@{$$self{chroms}}) {
			$self->spawn('check_vep',"$out/check_vep.$chr.done","$out/all.$chr.cons.vcf","$out/all.$chr.vcf.gz","$out/vep.$chr.done");
		}
	}
	$self->wait;
	# merge annotated VCFs (unfiltered)
    $self->set_limits( queue=>'small', memory=>1_000, runtime=>undef );
	foreach my $sample (sort keys %samples) {
		my $out = "$outdir/$sample";
		$self->spawn('concat_vcfs',"$out/all.cons.merged.vcf.gz",$out);
	}
	$self->wait;

	
	# option 1: do no filtering, VCF only
	if ($$self{output} eq 'vcf' && $$self{conseq_list} eq 'all' && !$$self{biotype}) {
		exit;
	}
	# option 2: no biotype filter, no conseq filter, 'both' or 'table'
	elsif (!$$self{biotype} && $$self{conseq_list} ne 'all' && $$self{output} ne 'vcf' ) {
		foreach my $sample (sort keys %samples) {
			my $out = "$outdir/$sample";
			$self->spawn('reformat_vcf',"$out/reformat.done","$out/all.cons.merged.vcf.gz",$$self{output},$out,$$self{add_DP4T},$$self{conseq_list});
		}
		$self->wait;
	}
	# option 3: filter by consequence/biotype, VCF only, both, or table only
	else {
		foreach my $sample (sort keys %samples) {
			my $out = "$outdir/$sample";
			$self->spawn('reformat_vcf',"$out/reformat.done","$out/all.cons.merged.vcf.gz",$$self{output},$out,$$self{add_DP4T},$$self{conseq_list},$$self{biotype_filter});
		}
		$self->wait;
	}
	close F;

	#$self->wait;
	# for use with 2 and 4 only:
	# optional: summary tab table with all tumour sample somatic calls (must run 2 or 3 first)
	if (($$self{output} eq 'both' || $$self{output} eq 'table') && $$self{summary_table} eq 'y') {
		my $list = "$$self{outdir}/mutect_table.list";
		#my @list;
		#foreach my $sample (sort keys %samples) {
		#	push @list, "$$self{outdir}/$sample/$sample-MuTect.txt";
		#}
		#open L, ">$list" or die;
		#	print L join "\n", @list;
		#close L;
		$self->spawn('mutectlist',$list,@mutect_list);
		$self->wait;
		$self->set_limits( queue=>'normal', memory=>3_000, runtime=>undef );
		$self->spawn('summary_table_sites',"$$self{outdir}/all_sites.annot",$$self{add_DP4T},$list);
		$self->wait;
		$self->set_limits( queue=>'small', memory=>1_000, runtime=>undef );
        for my $chr (@{$$self{chroms}}) {
			$self->spawn('summary_table',"$$self{outdir}/mutect_summary.$chr.txt",$chr,"$$self{outdir}/all_sites.annot",$list);
		}
		$self->wait;
		$self->spawn('concat_summary',"$$self{outdir}/summary_table-MuTect.txt",$$self{outdir});
		$self->wait;
	}
	# not implemented yet
    #$self->clean($outdir) unless (!$$self{do_clean});
	$self->cmd("touch $outdir/all_done");
    $self->all_done;

}

# clean not implemented yet
sub clean {
    my ($self,$outdir) = @_;
    $self->SUPER::clean($outdir);
    my $chunks = $self->get_chunks;
    for my $pop (keys %{$$self{pops}}) {
        for my $chunk (@$chunks) {
            my $chr  = $$chunk{chr};
            my $from = $$chunk{from};
            my $to   = $$chunk{to};
            for my $suffix (qw(samples vcf.gz vcf.gz.tbi)) {
                my $file = "$outdir/$pop/$chr/$chr:$from-$to.$suffix";
                unlink($file) unless !-e $file;
            }
        }
        for my $chr (@{$$self{chroms}}) {
            unlink("$outdir/$pop/$chr/concat.list") unless (!-e "$outdir/$pop/$chr/concat.list");
            rmdir("$outdir/$pop/$chr") unless ($$self{keep_bcfs});
        }
    }
}

sub save_config {
    my ($self) = @_;
    my $src = $$self{_config};
    my $dst = "$$self{outdir}/mutect_and_annotate.conf";
    if ( -e $dst && (stat($src))[9] <= (stat($dst))[9] ) { 
		return; 
	}
    if ( !-d $$self{outdir} ) { 
		$self->cmd("mkdir -p $$self{outdir}"); 
	}
    open(my $fh,'>',$dst) or $self->throw("$dst: $!");
    my $about = $$self{_about};
    $about =~ s/\n/\n# /g;
    print $fh "# $about";
    close($fh);
    $self->cmd("cat $src >> $dst");
}

sub cmd {
    my ($self,$cmd) = @_;
    $cmd =~ s/\n/ /g;
    return Utils::CMD($cmd,{verbose=>1});
}

sub symlinktarg {
	my ($self,$vcfout,$dir) = @_;
	my $cmd = "cd $$self{outdir}/$dir && ln -svf $dir.filt.targets.vcf.gz all.final.vcf.gz.tmp && ln -svf $dir.filt.targets.vcf.gz.tbi all.final.vcf.gz.tmp.tbi && cd $$self{outdir}";
    $self->cmd($cmd);
    rename("$vcfout.tmp.tbi","$vcfout.tbi") or $self->throw("rename $vcfout.tmp.tbi $vcfout.tbi: $!");
    rename("$vcfout.tmp",$vcfout) or $self->throw("rename $vcfout.tmp $vcfout: $!");

}

sub symlink {
	my ($self,$vcfout,$dir) = @_;
	my $cmd = "cd $$self{outdir}/$dir && ln -svf $dir.filt.vcf.gz all.final.vcf.gz.tmp && ln -svf $dir.filt.vcf.gz.tbi all.final.vcf.gz.tmp.tbi";
	$self->cmd($cmd);
	rename("$vcfout.tmp.tbi","$vcfout.tbi") or $self->throw("rename $vcfout.tmp.tbi $vcfout.tbi: $!"); 
	rename("$vcfout.tmp",$vcfout) or $self->throw("rename $vcfout.tmp $vcfout: $!"); 
}

sub mutectlist {
	my ($self,$list,@samples) = @_;
	my @list;
	foreach my $sample (@samples) {
		my @filename = `dir $$self{outdir}/$sample/*-MuTect.txt`;
		if (!@filename || (@filename && @filename > 1)) {
			die "Problem with *-MuTect.txt file in  $$self{outdir}/$sample\n";
		}
		chomp @filename;
		push @list, $filename[0];
		#push @list, "$$self{outdir}/$sample/$sample-MuTect.txt";
	}
	open L, ">$list" or die;
		print L join "\n", @list;
	close L;

}	

sub mutect_call {
    my ($self,$outfile,$prefix,$tbam,$nbam,$chr) = @_;
#java -Xmx2g -jar /software/team113/algorithms/muTect/muTect-1.1.4.jar --analysis_type MuTect --reference_sequence /lustre/scratch105/vrpipe/refs/human/ncbi37/human_g1k_v37.fasta --dbsnp  /lustre/scratch105/vrpipe/refs/human/ncbi37/dbsnp/dbsnp_141.b37.vcf.gz --input_file:normal ./BAMS/A04_LCL/A04_LCL.bam --input_file:tumor ./BAMS/A04/A04.bam --out test.stats --coverage_file  test.cov   --cosmic /lustre/scratch107/user/kw10/mouse/SEQCAP_Identification_of_drug_resistance_genes_in_melanoma/BAMS_REMAPPED/MUTECT/CosmicCodingMuts.vcf.gz --test.vcf
	my $vcf = $chr ? "$prefix.$chr.vcf" : "$prefix.vcf";
	my $out = $chr ? "$prefix.mutect_stats.$chr.txt" : "$prefix.mutect_stats.txt";
	my $wig= $chr ? "$prefix.coverage.$chr.wig" : "$prefix.coverage.wig";

	my $cmd = "$$self{java} -Xmx2g -jar $$self{mutect} --enable_extended_output --analysis_type MuTect --reference_sequence $$self{reference} --vcf $vcf --input_file:normal $nbam  --input_file:tumor $tbam --out $out --coverage_file $wig ";
	if ( $$self{snpfile} ) {
		$cmd .= "--dbsnp $$self{snpfile} ";
	}
	if ( $$self{cosmic} )  {
		$cmd .= "--cosmic $$self{cosmic} ";
	}
	if ( $chr ) { 
		$cmd .= "-L $chr ";
	}
	if ( $$self{mutect_opt} ) {
		$cmd .= "$$self{mutect_opt} ";
	}
	$self->cmd($cmd);
	$self->cmd("touch $outfile");
}

# sort and concat split mutect output
# stats file and vcf file only
# remove wig files and idx

sub mutect_merge_vcf {
	my ($self,$vcfout,$outdir,$prefix) = @_;
	my $chrlist = join (" ", @{$$self{chroms}});
	print STDOUT "CHR: $chrlist\n";
	my $getlist = "for f in $chrlist; do echo $outdir/$prefix.\$f.vcf; done > $outdir/vep.mergelist";
	system($getlist) == 0 || die "$getlist failed : $!";
	my $cmd = "$$self{bcftools} concat -f $outdir/vep.mergelist -O z > $vcfout.tmp && tabix -p vcf $vcfout.tmp";
	$self->cmd($cmd);
	rename("$vcfout.tmp.tbi","$vcfout.tbi") or $self->throw("rename $vcfout.tmp.tbi $vcfout.tbi: $!"); 
	rename("$vcfout.tmp",$vcfout) or $self->throw("rename $vcfout.tmp $vcfout: $!"); 
}

sub mutect_merge_stats {
	my ($self,$statsout,$outdir,$prefix) = @_;
	my $chrlist = join (" ", @{$$self{chroms}});
	print STDOUT "CHR: $chrlist\n";
	my $getheader = "head -n 2 $outdir/$prefix.mutect_stats.$$self{chroms}[0].txt > $statsout.tmp";
	system($getheader) == 0 || die "$getheader failed  $!";
	my $getlist = "for f in $chrlist; do echo $outdir/$prefix.mutect_stats.\$f.txt; done > $outdir/mutect.mergelist";
	system($getlist) == 0 || die "$getlist failed : $!";
	my $cmd = qq(xargs cat < $outdir/mutect.mergelist  | grep -v ^# | perl -ne 'print unless /^contig\\tposition/' | sort -uV >> $statsout.tmp);
	$self->cmd($cmd);
	##rename("$statsout.tmp.tbi","$statsout.tbi") or $self->throw("rename $statsout.tmp.tbi $statsout.tbi: $!"); 
	rename("$statsout.tmp",$statsout) or $self->throw("rename $statsout.tmp $statsout: $!"); 
}

sub add_dp4 {
	my ($self,$vcfout,$remove_reject,$statfile,$vcfin,$dir,$sample) = @_;
	my $prefix = "$dir/$sample";
	my $header = '"##INFO=<ID=DP4T,Number=4,Type=Integer,Description=\"Number of tumour REF for, REF rev, ALT for, ALT ref bases\">"';
	print STDOUT "PATH: echo $header $prefix/dp4header";
	`echo $header > $prefix/dp4header`;
	my @statsheader = `head -n 2 $statfile | tail -n 1`;
	my ($contig,$pos,$strandcounts);
	my $i=1;
	foreach (split /\t/, $statsheader[0]) {
		$contig = $i if $_ eq 'contig';
		$pos = $i if $_ eq 'position';
		$strandcounts = $i if $_ eq 'strand_bias_counts';
		$i++;
	}
	die if !$contig || !$pos || !$strandcounts;
	my $cmd;
	if ($remove_reject eq 'n') {
		# different format for version 1.1.7
		$cmd = "cut -f $contig,$pos,$strandcounts $statfile | grep -v ^contig | grep -v ^#  | sed 's/(//; s/)//' > $prefix/DP4T.annot.tab && bgzip -f $prefix/DP4T.annot.tab && tabix -f -s 1 -b 2 -e 2 $prefix/DP4T.annot.tab.gz && bcftools annotate -a $prefix/DP4T.annot.tab.gz -c CHROM,POS,DP4T  -h $prefix/dp4header $vcfin | bgzip -fc > $vcfout.tmp && tabix -f -p vcf $vcfout.tmp && rm -f $prefix/dp4header $prefix/DP4T.annot.tab*";
	}
	else {
		# different format for version 1.1.7
		$cmd = "cut -f $contig,$pos,$strandcounts $statfile | grep -v ^contig | grep -v ^#  | sed 's/(//; s/)//' > $prefix/DP4T.annot.tab && bgzip -f $prefix/DP4T.annot.tab && tabix -f -s 1 -b 2 -e 2 $prefix/DP4T.annot.tab.gz && bcftools annotate -a $prefix/DP4T.annot.tab.gz -c CHROM,POS,DP4T  -h $prefix/dp4header $vcfin |  awk '/^#/ || \$7!~/REJECT/' | bgzip -fc > $vcfout.tmp && tabix -f -p vcf $vcfout.tmp && rm -f $prefix/dp4header $prefix/DP4T.annot.tab*";
	}
	$self->cmd($cmd);
	rename("$vcfout.tmp.tbi","$vcfout.tbi") or $self->throw("rename $vcfout.tmp.tbi $vcfout.tbi: $!"); 
	rename("$vcfout.tmp",$vcfout) or $self->throw("rename $vcfout.tmp $vcfout: $!"); 
	
}

sub filter_reject {
	my ($self,$vcfout,$vcfin) = @_;
	my $out = $$self{outdir};
	my $cmd = "awk '/^#/ || \$7!~/REJECT/' $vcfin | bgzip -fc > $vcfout.tmp && tabix -f -p vcf $vcfout.tmp";
	$self->cmd($cmd);
	rename("$vcfout.tmp.tbi","$vcfout.tbi") or $self->throw("rename $vcfout.tmp.tbi $vcfout.tbi: $!"); 
	rename("$vcfout.tmp",$vcfout) or $self->throw("rename $vcfout.tmp $vcfout: $!"); 
}

sub filter_targets {
	my ($self,$vcfout,$vcfin) = @_;
	my $cmd = "$$self{bedtools} intersect -a $vcfin -b $$self{targets} -header | bgzip -c > $vcfout.tmp && tabix -f -p vcf $vcfout.tmp ";
	$self->cmd($cmd);
	rename("$vcfout.tmp.tbi","$vcfout.tbi") or $self->throw("rename $vcfout.tmp.tbi $vcfout.tbi: $!"); 
	rename("$vcfout.tmp",$vcfout) or $self->throw("rename $vcfout.tmp $vcfout: $!"); 
}

sub split_vcf {
	my ($self,$vcfout,$vcfin,$chr) = @_;
	my $cmd = "tabix -h $vcfin $chr | bgzip -c > $vcfout.tmp && tabix -p vcf $vcfout.tmp"; 
	$self->cmd($cmd);
	rename("$vcfout.tmp.tbi","$vcfout.tbi") or $self->throw("rename $vcfout.tmp.tbi $vcfout.tbi: $!"); 
	rename("$vcfout.tmp",$vcfout) or $self->throw("rename $vcfout.tmp $vcfout: $!"); 
}
sub run_vep {
	my ($self,$outfile,$vcfout,$vcfin) = @_;
	# run vep by chrom
	# vep
	# vep outputs an unzipped vcf txt file
	my $cmd = "perl -I $$self{ensembl_api} $$self{vep} -i $vcfin --db_version $$self{ens_version}  -t SO --format vcf -o $vcfout.tmp --force_overwrite --cache --dir $$self{vep_cache} --buffer 20000 --species $$self{species} --offline --symbol --biotype --vcf --sift s ";
	$cmd .= "--no_stats " if $$self{vep_stats} eq 'no' || $$self{vep_stats} eq 'n' ;
	$cmd .= "--assembly $$self{assembly} " if $$self{assembly};
	$cmd .= "--hgvs --shift_hgvs 1 " if $$self{add_hgvs};
	if ( $$self{vep_fasta} ) { $cmd .= " --fasta $$self{reference}"; }
	$self->cmd($cmd);
	if ( -e "$vcfout.tmp") {
		rename("$vcfout.tmp",$vcfout) or $self->throw("rename $vcfout.tmp $vcfout: $!"); 
	}
	$self->cmd("touch $outfile");
}

sub check_vep {
	my ($self,$outfile,$infile,$oldfile,$donefile) = @_;
	if (! -e $infile) {
		my $cmd = "zcat $oldfile > $infile";
		system($cmd) == 0 || die "$cmd failed : $!";
		system("touch $donefile") == 0 || die "touch $donefile failed : $!";
	}
	my $cmd2 = "bgzip $infile && tabix -p vcf $infile.gz";
	$self->cmd($cmd2);
	$self->cmd("touch $outfile");

#	rename("$outfile.tmp.gz.tbi","$outfile.tbi") or $self->throw("rename $outfile.tmp.gz.tbi $outfile.tbi: $!"); 
#	rename("$outfile.tmp.gz",$outfile) or $self->throw("rename $outfile.tmp.gz $outfile: $!"); 

}

# sort and concat annotated VCFs
sub concat_vcfs {
	my ($self,$vcfout,$outdir) = @_;
	my $chrlist = join (" ", @{$$self{chroms}});
	print STDOUT "CHR: $chrlist\n";
	my $getlist = "for f in $chrlist; do echo $outdir/all.\$f.cons.vcf.gz; done > $outdir/vep.mergelist";
	system($getlist) == 0 || die "$getlist failed : $!";
	my $cmd = "$$self{bcftools} concat -f $outdir/vep.mergelist -O z > $vcfout.tmp && tabix -p vcf $vcfout.tmp";
	$self->cmd($cmd);
	rename("$vcfout.tmp.tbi","$vcfout.tbi") or $self->throw("rename $vcfout.tmp.tbi $vcfout.tbi: $!"); 
	rename("$vcfout.tmp",$vcfout) or $self->throw("rename $vcfout.tmp $vcfout: $!"); 
}

sub reformat_vcf {
	my ($self,$outfile,$vcfin,$outputformats,$outdir,$which,$conseq,$biotype) = @_;
	$which = $which eq 'y' ? 'extended' : 'default';
	my $cmd = "zcat $vcfin | $$self{reformat_vcf} -f $conseq -m $which -o $outputformats -d $outdir ";
	if ($biotype) {
		$cmd .= "-b $biotype";
	}
	$self->cmd($cmd);
	$self->cmd("touch $outfile");
}

sub summary_table_sites {
	my ($self,$outfile,$dp4,$list) = @_;
	my $header = "#CHROM\tPOS\tID\tREF\tALT\tENS_ID\tGENE_SYMBOL\tCONSEQUENCE\tCDS_POS\tPROT_POS\tAA_CHANGE\tTRANSCRIPTS\tTRANS_BIOTYPE\tSIFT";
	my $cmd;
	if ($dp4 eq 'n') {
		#$cmd = "echo -e \"$header\" > $$self{outdir}/all_sites.annot.tmp && for f in `cat $list`; do cut -f 1-5,11- \$f | grep -v ^# ; done |  sort -uV >> $$self{outdir}/all_sites.annot.tmp";
		$cmd = "echo -e \"$header\" > $$self{outdir}/all_sites.annot.tmp && xargs cat < $list | grep -v ^# | cut -f 1-5,11- | sort -uV >> $$self{outdir}/all_sites.annot.tmp";
	}
	else {
		#$cmd = "echo -e \"$header\" > $$self{outdir}/all_sites.annot.tmp && for f in `cat $list`; do cut -f 1-5,12- \$f | grep -v ^# ; done |  sort -uV >> $$self{outdir}/all_sites.annot.tmp";
		$cmd = "echo -e \"$header\" > $$self{outdir}/all_sites.annot.tmp && xargs cat < $list | grep -v ^# | cut -f 1-5,12- | sort -uV >> $$self{outdir}/all_sites.annot.tmp";
	}
	$self->cmd($cmd);
	rename("$$self{outdir}/all_sites.annot.tmp",$outfile) or $self->throw("rename $$self{outdir}/all_sites.annot.tmp $outfile: $!"); 
}

sub summary_table {
	my ($self,$outfile,$chr,$sitelist,$filelist) = @_;
	my $cmd2 = "awk '\$1~/^#/ || \$1==\"$chr\"' $sitelist | $$self{mutect_summary} $filelist > $$self{outdir}/mutect_summary.$chr.tmp";
	$self->cmd($cmd2);
	rename("$$self{outdir}/mutect_summary.$chr.tmp",$outfile) or $self->throw("rename $$self{outdir}/mutect_summary.$chr.tmp $outfile: $!"); 
}

sub concat_summary {
	my ($self,$outfile,$out) = @_;
	my @list2;
	for my $chr (@{$$self{chroms}}) {
		push @list2, "$$self{outdir}/mutect_summary.$chr.txt";
	}
	my $list2 = "$$self{outdir}/mutect_summary.list";
	open L2, ">$list2" or die;
		print L2 join "\n", @list2;
	close L2;
	my @list = `cat $list2`;
	chomp @list;
	my $catlist = join (" ", @list);
	my $cmd = "grep ^# $list[0] > $out/mutect_summary.tmp  && cat $catlist | grep -v ^# | sort -V >> $out/mutect_summary.tmp";
	#my $cmd = "grep ^# $out/mutect_summary.$chrom[0].txt > $out/mutect_summary.tmp && for f in $chrlist; do cat $out/mutect_summary.\$f.txt | grep -v ^# >> $out/mutect_summary.tmp ; done";
	$self->cmd($cmd);
	rename("$out/mutect_summary.tmp",$outfile) or $self->throw("rename $out/mutect_summary.tmp $outfile: $!"); 
}

#sub filter_vcf {
#	my ($self,$vcf,$cons,$biotype) = @_;
#	my $cmd;
#	if ($cons && $biotype) {
#		$cmd = "zcat $vcf | awk '/^#|$cons/' | awk '/^#|$biotype/' |bgzip -c > $vcfout && tabix -p vcf $vcfout.gz";
#		$self->cmd($cmd); 
#
#	}
#	elsif ($cons || $biotype) {
#		my $filt = $cons ? $cons : $biotype;
#		$cmd = "zcat $vcf | awk '/^#|$filt/' | bgzip -c > $vcfout && tabix -p vcf $vcfout.gz";
#		$self->cmd($cmd); 
#	}
#	
#
#}
