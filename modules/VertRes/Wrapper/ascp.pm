
=head1 NAME

VertRes::Wrapper::ascp - wrapper for fastqcheck

=head1 SYNOPSIS

use VertRes::Wrapper::ascp;

my $wrapper = VertRes::Wrapper::ascp->new();

$wrapper->upload_file('file.local.bam', 'file.remote.bam');
$wrapper->upload_directory('/local/dir/path/', '/remote/dir/path/');
$wrapper->upload_fofn('local.fofn', 'remote/path/');

my $wrapper = VertRes::Wrapper::ascp->new(dropbox => 'g1k'); # for g1k dropbox (default)
my $wrapper = VertRes::Wrapper::ascp->new(dropbox => 'g1k'); # for g1k dropbox (default)

Otherwise use the following to specify the connection details for the dropbox
my $wrapper = VertRes::Wrapper::ascp->new(username => 'username', server => 'server', password => 'password');


=head1 DESCRIPTION

Transfer files to or from an Aspera dropbox.

=head1 AUTHOR

Shane McCarthy: sm15@sanger.ac.uk

=cut

package VertRes::Wrapper::ascp;      

use strict;
use warnings;
use VertRes::IO;
use Net::SSH::Expect;
use File::Basename;
use File::Listing;
use Cwd;

use base qw(VertRes::Wrapper::WrapperI);

my $server = 'fasp.sra.ebi.ac.uk';
my $username = 'g1k-drop-si';
my $password = 'Rn3znNB5';


=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::ascp->new();
 Function: Create a VertRes::Wrapper::ascp object.
 Returns : VertRes::Wrapper::ascp object
 Args    : quiet => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args,
    								exe => '~sm15/.aspera/connect/bin/ascp',
    								run_method => 'write_to',
									params => [qw(l m M k K)],
									switches => [qw(A Q T d p q v)]);
    
    return $self;
}

=head2 version

 Title   : version
 Usage   : my $version = $obj->version();
 Function: Returns the program version.
 Returns : Aspera scp version or 0 if not found.
 Args    : n/a

=cut

sub version {
    my $self = shift;
    
    my $exe = $self->exe;
    open(my $fh, "$exe -A 2>&1 |") || $self->throw("Could not start $exe");
    my $version = 0;
    while (<$fh>) {
        if (/version (\S+)/) {
            $version = $1;
            last;
        }
    }
    close($fh);
    
    return $version;
}


=head2 upload_file

 Title   : upload_file
 Usage   : $ascp->upload_file($source_file, $dest_file);
 Function: Upload a file. A specified number (default 2) of attempts will be made. 
 		   Upon success the file .$source_file.aspera will be created so that we know 
 		   this file has been uploaded and we will not upload again if it happens
 		   to be removed from the dropbox.
 Returns : boolean
 Args    : paths to input and output files

=cut


sub upload_file {
	my ($self, $source_file, $dest_file, %opts) = @_;

	my $in_dropbox = $self->match_filesize($source_file, $dest_file);
	
	my $log_file = ".$source_file.aspera";
	if (defined $opts{refresh}) {
		unlink $log_file if (-e $log_file);
	}
	
	if ($in_dropbox || -e $log_file) {
		print "Already uploaded... skipping.\n";
		#if $opts{verbose}
		return 1;
	}
	
	my $attempt = 0;
	my $max_attempt = 2;
	if (defined $opts{max_attempt}) {
		$max_attempt = $opts{max_attempt};
	}
	
	while (!$in_dropbox && $attempt < $max_attempt) {
		$in_dropbox = $self->_upload($source_file, $dest_file);
		++$attempt;
	}
	
	if ($in_dropbox) {
		print "Upload of $source_file successful.\n";
		#if $opts{verbose}
# 		system("touch $log_file \n");
		return 1;
	} else {
		print "Upload of $source_file unsuccessful.\n";
		#if $opts{verbose}
		return 0;
	}
}

=head2 _upload

 Title   : upload_file
 Usage   : $obj->upload_file($input_file, $output_file);
 Function: Upload a file.
 Returns : boolean
 Args    : paths to input and output files

=cut

sub _upload {
	my ($self, $source_file, $dest_file) = @_;
	
    $self->_set_params_and_switches_from_args(k => 1, Q => 1, T => 1, m => '50M', l => '400M');
		
	my $pipe = $self->run($source_file, "$username\@$server:$dest_file");
	print $pipe "$password\n";
	close $pipe;
	
	return $self->match_filesize($source_file, $dest_file);
}

=head2 upload_directory

 Title   : upload_directory
 Usage   : $obj->upload_directory($source_dir, $out_dir);
 Function: Upload everything in a directory.
 Returns : boolean
 Args    : paths to input and output files

=cut

sub upload_directory {
	my ($self, $source_dir, $dest_dir) = @_;
		
	opendir(my $dh, $source_dir) || die "can't opendir $source_dir: $!";
    my @source_files = readdir($source_dir);
    closedir $dh;

	my $number_of_files = scalar @source_files;
	my $number_uploaded = 0;
	
	foreach my $source_file (@source_files) {
		my $dest_file = $dest_dir.basename($_);
		$number_uploaded++ if $self->upload_file($source_file, $dest_file);
	}
	
	my $match = $number_of_files == $number_uploaded ? 1 : 0;
	
	return $match;
}

=head2 upload_fofn

 Title   : upload_fofn
 Usage   : $obj->upload_fofn($input_file, $output_file);
 Function: Upload a all files in a file of filenames.
 Returns : boolean. True (1) if all files in the fofn have uploaded.
 Args    : paths to input and output files

=cut

sub upload_fofn {
	my ($self, $fofn, $dest_dir) = @_;
	
	open my $fh, "<$fofn";
	$fh || $self->throw("failed to open filehandle for fofn '$fofn'");
	
	my $number_of_files = 0;
	my $number_uploaded = 0;
	
	while( <$fh> ) {
		chomp;
		$number_of_files++;
		my $dest_file = $dest_dir.basename($_);
		$number_uploaded++ if $self->upload_file($_, $dest_file);
	}
	close $fh;
	
	my $match = $number_of_files == $number_uploaded ? 1 : 0;
	
	return $match;
}




=head2 match_filesize

 Title   : match_filesize
 Usage   : my $version = $obj->version();
 Function: Compares filesizes of source and destination by 
 	logging in to aspera dropbox via SSH and ls-ing the file size.
 Returns : Boolean - true (1) if filesizes are the same.
 Args    : source file, destination file

=cut


sub match_filesize {
	my ($self, $source, $target) = @_;
	
    my $ssh = Net::SSH::Expect->new (
        host => $server, 
        password=> $password, 
        user => $username, 
        raw_pty => 1
    );
    my $login_output = $ssh->login();
	if ($login_output !~ /aspsh>/) {
		die "Login has failed. Login output was $login_output";
	}
	
	my $target_dir = File::Basename::dirname($target);
    my $out = $ssh->exec("ls -l $target_dir/");
    $ssh->close();
	
	my $orig_size = -s $source;
	
	my $new_size = 0;
	for (File::Listing::parse_dir("$out\n")) {
		my ($name, $type, $size, $mtime, $mode) = @$_;
		if ($name eq File::Basename::basename($target)) {
			$new_size = $size;
			last;
		}
	}

	my $match = $orig_size == $new_size ? 1 : 0;
	
	return $match;
}




