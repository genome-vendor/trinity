#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib/");
use SAM_reader;
use SAM_entry;
use Fasta_reader;
use EM;
use Getopt::Long qw(:config no_ignore_case bundling);
use Data::Dumper;

## written by bhaas, equations for EM based on work by Moran Yassour and Nir Friedman

my $usage = <<__EOUSAGE__;

#####################################################################################################################
#
#  ## Required:
#
#  --sam      *read name* sorted sam file
#  --fasta    transcripts.fasta
#
#  ## Optional
#
#  --max_iterations   maximal number of EM iterations (default: 1000)
#  --min_delta_ML     EM stops when difference between subsequent max likelihood computations is below this value. (default: 0.001)
#
#####################################################################################################################

__EOUSAGE__

	;


my ($sam_file, $transcripts_fa);

my $MAX_ITERATIONS = 1000;
my $MIN_DELTA_ML = 0.001;


&GetOptions( 'sam=s' => \$sam_file,
			 'fasta=s' => \$transcripts_fa,
			 'max_iterations=i' => \$MAX_ITERATIONS,
			 'min_delta_ML=f' => \$MIN_DELTA_ML,
			 );

unless ($sam_file && $transcripts_fa) {
	die $usage;
}


unless ($sam_file =~ /\.sam$/) {
	die "Error, sam file requires .sam extension";
}


main: {

	
	print STDERR "-processing FASTA sequences: $transcripts_fa\n";
	
	my $fasta_reader = new Fasta_reader($transcripts_fa);
		
	my %transcript_seqs = $fasta_reader->retrieve_all_seqs_hash();

	## examine each component's data at a time.
	
	my $sam_reader = new SAM_reader($sam_file);
	
	my $em_obj = EM->new(\%transcript_seqs);	

	
	my %transcripts;
	my $curr_read_name = "";

	print STDERR "-processing SAM file: $sam_file\n";

	my $read_counter = 0;
	while (my $sam_entry = $sam_reader->get_next()) {
				
		my $trans_id = $sam_entry->get_scaffold_name();
		my $read_name = $sam_entry->get_read_name();
		
		if ($read_name ne $curr_read_name) {
			$read_counter++;
			print STDERR "\r[$read_counter] $read_name     " if $read_counter % 1000 == 0;
			
			if (%transcripts) {
				$em_obj->add_read(keys %transcripts);
				%transcripts = ();
			}
		}
		$curr_read_name = $read_name;
		$transcripts{$trans_id} = 1;
		
	}
	
	# get last ones
	if (%transcripts) {
		$em_obj->add_read(keys %transcripts);
	}
		
	print STDERR "\n\n-running EM\n\n";
	
	$em_obj->run($MAX_ITERATIONS, $MIN_DELTA_ML);

	$em_obj->print_results();
	

	exit(0);
	
}

