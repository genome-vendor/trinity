#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use Getopt::Long qw(:config no_ignore_case bundling);

my $usage = <<_EOUSAGE_;

#################################################################################
#
#  --reads <string>  filename.fasta          
#
#  --minK <int>      Default: 19   (must be odd number)
#  --maxK <int>      Default: 29   (must be odd number, avoid palandromes)
#
#  --CPU             Default: 1
#
#  --iworm_opts <string>   additional options to give to inchworm, other than -K!
#
###################################################################################


_EOUSAGE_

	;


my $reads_file;
my $minK = 19;
my $maxK = 29;
my $CPU = 1;
my $iworm_opts = "";

&GetOptions (
	"reads=s" => \$reads_file,
	"minK=i" => \$minK,
	"maxK=i" => \$maxK,
	"CPU=i" => \$CPU,
	"iworm_opts=s" => \$iworm_opts,
	);


unless ($reads_file) {
	die $usage;
}

if ($maxK < $minK) {
	die $usage;
}

if ($maxK % 2 == 0 || $minK % 2 == 0) {
	die $usage;
}

if ($iworm_opts =~ /-K /) {
	die "Error, cannot specify -K in iworm_opts";
}



my $log_dir = "iworm_log.$$";
mkdir($log_dir) or die "Error, cannot mkdir $log_dir";

## This is very cheap 'parallel-computing' !!!  :)

main: {

	my %job_tracker;
	
	my $num_running = 0;
	
	for (my $i = $minK; $i <= $maxK; $i += 2) {
		
		if ($num_running < $CPU) {
			
			$num_running++;
			my $child = fork();

			if ($child) {
				# parent
				$job_tracker{$child} = $i;
			}
			else {
				# child:
				my $ret = &run_iworm($i);
				exit($ret);
			}
		}
		
		if ($num_running >= $CPU) {
			wait();
			$num_running--;
		}
	}

	while (wait() != -1) { };

	if (&all_succeeded(\%job_tracker)) {

		&run_iworm_reassembly(\%job_tracker);
		
		`rm -rf $log_dir`;
		
	}


	exit(0);
}

	
####
sub run_iworm {
	my ($kmer) = @_;

	my $cmd = "$FindBin::Bin/inchworm --reads $reads_file --run_inchworm -K $kmer $iworm_opts > $reads_file.K$kmer.iworm";

	print "\nRUNNING: $cmd\n\n";
	
	my $ret = 0;

	## don't rerun it if result file already exists.
	if (-s "$reads_file.K$kmer.iworm") {
		print "** Warning: output file: $reads_file.K$kmer.iworm  already exists.... Keeping it rather than rerunning inchworm.\n";
	}
	else {
		$ret = system($cmd);
	}
	
	if ($ret) {
		print STDERR "Error, command: $cmd died with ret $ret";
	}

	open (my $log_fh, ">$log_dir/$kmer.ret") or die "Error, cannot write to log file for $kmer.ret";
	print $log_fh $ret;
	close $log_fh;


	return($ret);
}


####
sub all_succeeded {
	my ($job_tracker_href) = @_;

	my @kmers = values %$job_tracker_href;

	my $OK = 1;

	foreach my $kmer (@kmers) {
		my $log_file = "$log_dir/$kmer.ret";
		my $ret_val = `cat $log_file`;
		chomp $ret_val;

		unless (-s $log_file && $ret_val == 0) {
			$OK = 0;
			print STDERR "Error, kmer run $kmer failed.\n";
		}
	}

	return($OK);
}


####
sub run_iworm_reassembly {
	my ($job_tracker_href) = @_;

	my @kmers = sort {$a<=>$b} values (%$job_tracker_href);

	my $cmd = "cat ";
	foreach my $kmer (@kmers) {
		$cmd .= "$reads_file.K$kmer.iworm ";
	}

	$cmd .= " > all_indiv_iworms.fasta";

	my $ret = system($cmd);


	if ($ret) {
		die "Error, command $cmd died with ret $ret";
	}


	## run iworm now in reassembly mode:

	my $min_kmer = shift @kmers;
	my $max_kmer = pop @kmers;
	
	$cmd = "$FindBin::Bin/inchworm --reads all_indiv_iworms.fasta --run_inchworm --reassembleIworm $iworm_opts > varyK.$min_kmer-$max_kmer.reassembly.iworm";

	$ret = system($cmd);

	if ($ret) {
		die "Error, command $cmd died with ret $ret";
	}

	else {
		return;
	}
}


		
