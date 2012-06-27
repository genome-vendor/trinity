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

###########################################################################################################################################
#
#  ## Required:
#
#  --coord_sorted_sam      coord-sorted sam  (  sort -k3,3 -k4,4n file.sam > coordSorted.sam )
#  --name_sorted_sam      read name-sorted sam  (   sort -k1,1 -k3,3 file.sam > nameSorted.sam )
#  --fasta    Trinity.fasta
#
#  ## Optional
#
#  --out_prefix           prefix for output files. (default: Trinity.refined, resulting in Trinity.refined.fasta and Trinity.refined.expression)
#
#  --retain_min_percent_component_expressed (default: 1)
#
#  Phase 1 (component pre-screen):
#  --max_iterations_phase1    component-based pre-screen, maximal number of EM iterations (default: 1000)
#  --no_require_unique        by default, requires each transcript in each component to have at least one unique read supporting it.
#  --aggressive               rapidly removes transcripts that have no unique read mapping.  By default, deletes one at a time and reexamines uniqueness after each deletion.
#
#  Phase 2: (full data set simultaneously)
#  --max_iterations_phase2    final estimation of abundance, max number of EM interations (default: 1000)
#
#  Both phases>
#  --min_delta_ML     EM stops when difference between subsequent max likelihood computations is below this value. (default: 0.01)
#
############################################################################################################################################


__EOUSAGE__

	;


my ($coord_sorted_sam_file, $name_sorted_sam_file, $transcripts_fa);

my $MAX_ITERATIONS_PHASE1 = 1000;
my $MAX_ITERATIONS_PHASE2 = 1000;
my $MIN_DELTA_ML = 0.01;
my $RETAIN_MIN_PERCENT_COMPONENT_EXPRESSED = 1.0;
my $help_flag;
my $OUT_PREFIX = "Trinity.refined";
my $NO_REQUIRE_UNIQUE_READ_PHASE1 = 0;
my $AGGRESSIVE = 0;

&GetOptions( 'coord_sorted_sam=s' => \$coord_sorted_sam_file,
			 'name_sorted_sam=s' => \$name_sorted_sam_file,
			 'fasta=s' => \$transcripts_fa,
			 'max_iterations_phase1=i' => \$MAX_ITERATIONS_PHASE1,
			 'max_iterations_phase2=i' => \$MAX_ITERATIONS_PHASE2,
			 'min_delta_ML=f' => \$MIN_DELTA_ML,
			 'retain_min_percent_component_expressed=f' => \$RETAIN_MIN_PERCENT_COMPONENT_EXPRESSED,
			 'h' => \$help_flag,
			 'out_prefix=s' => \$OUT_PREFIX,
			 'no_require_unique' => \$NO_REQUIRE_UNIQUE_READ_PHASE1,
			 'aggressive' => \$AGGRESSIVE,
			 );



if ($help_flag) {
	die $usage;
}

unless ($coord_sorted_sam_file && $name_sorted_sam_file && $transcripts_fa) {
	die $usage;
}


main: {

	
	print STDERR "-processing FASTA sequences: $transcripts_fa\n";
	
	my $fasta_reader = new Fasta_reader($transcripts_fa);
		
	my %transcript_seqs = $fasta_reader->retrieve_all_seqs_hash();

	## examine each component's data at a time.
		
	## Phase 1: examine expression on per-component basis.

	print STDERR "\n\n\n#######   Phase 1, Per-Component Analysis/Filtering ##########\n\n";
	my %transcripts_passed = &examine_EM_per_component($coord_sorted_sam_file, \%transcript_seqs);  ## meet min % isoform expressed requirement in per-component analysis

	## Phase 2: report final transcripts and abundance estimates.
	print STDERR "\n\n\n########   Phase 2, Analysis Exploring All Reads and All Components Simutaneously  ###########\n\n\n";
	&compute_abundance_estimates_filtered_Trinity(\%transcripts_passed, $name_sorted_sam_file, \%transcript_seqs);

	print "\n\nFinished.  See output files: $OUT_PREFIX.\*\n\n";
	

	exit(0);
}


#### 
sub examine_EM_per_component {
	my ($coord_sorted_sam_file, $transcript_seqs_href) = @_;

	
	my %transcripts_OK;
	
	my $curr_component = undef;
	my %read_to_trans;

	print STDERR "-processing coord_sorted_sam_file: $coord_sorted_sam_file\n";
		
	my $sam_reader = new SAM_reader($coord_sorted_sam_file);
	
	my $read_counter = 0;
	while (my $sam_entry = $sam_reader->get_next()) {
		$read_counter++;
		if ($read_counter % 1000 == 0) {
			print STDERR "\r[$read_counter] sam entries read.   ";
		}
		
		my $trans_id = $sam_entry->get_scaffold_name();
		my $read_name = $sam_entry->get_read_name();
		
		$trans_id =~ /^(comp\d+)_/ or die "Error, cannot parse Trinity component number from $trans_id";
		my $component = $1;
		if ($curr_component && $curr_component ne $component) {
			if (%read_to_trans) {
				&run_component_EM(\%read_to_trans, \%transcripts_OK, $transcript_seqs_href, $curr_component);
				%read_to_trans = (); # reinit
			}
		}
		$curr_component = $component;
		$read_to_trans{$read_name}->{$trans_id} = 1;
	}

	if (%read_to_trans) {
		&run_component_EM(\%read_to_trans, \%transcripts_OK, $transcript_seqs_href, $curr_component); # get last ones
	}
		
	return(%transcripts_OK);
}


####
sub run_component_EM {
	my ($read_to_trans_href, $transcripts_OK_href, $transcript_seqs_href, $curr_component) = @_;
	

	my $em_obj = EM->new($transcript_seqs_href);


	unless ($NO_REQUIRE_UNIQUE_READ_PHASE1) {
		$read_to_trans_href = &require_unique_read_mapping($read_to_trans_href);
	}
	
	foreach my $read_name (keys %$read_to_trans_href) {
		
		my $transcripts_href = $read_to_trans_href->{$read_name};
		
		$em_obj->add_read(keys %$transcripts_href);
	
	}
	
	my @all_transcripts = $em_obj->_get_all_transcripts();
	
	if (scalar (@all_transcripts) == 1) {
		## no reason to run EM.
		$transcripts_OK_href->{ $all_transcripts[0] } = 1;
		return;
	}
	
	print STDERR "-running EM on " . scalar(@all_transcripts) . " Trinity transcripts of component: " . $curr_component . "\n";
		
	$em_obj->run(max_iterations => $MAX_ITERATIONS_PHASE1, 
				 min_delta_ML => $MIN_DELTA_ML);
	
	$em_obj->report_results();
	
	my @results = $em_obj->get_results();
   
	my $sum_fpkm = 0;
	foreach my $result (@results) {
		$sum_fpkm += $result->{FPKM};
	}
	
	foreach my $result (@results) {
		my $percent_expressed = $result->{FPKM} / $sum_fpkm * 100;
		my $transcript = $result->{trans_id};
		if ($percent_expressed >= $RETAIN_MIN_PERCENT_COMPONENT_EXPRESSED) {
			$transcripts_OK_href->{ $transcript } = 1;
			print "\t-retaining $transcript (" . sprintf("%.2f", $percent_expressed) . "% of component-level expression)\n";
		}
		else {
			print "\t-rejecting $transcript (" . sprintf("%.2f", $percent_expressed) . "% of component-level expression)\n";
			
		}
	}
	
	return;
	
}

####
sub compute_abundance_estimates_filtered_Trinity {
	my ($transcripts_passed_href, $name_sorted_sam_file, $transcript_seqs_href) = @_;

	print STDERR "-processing name_sorted_sam_file: $name_sorted_sam_file\n\n";
   
    
	my $em_obj;

	my $filtered_transcripts_flag = 1;
	while ($filtered_transcripts_flag) {
		
		my $sam_reader = new SAM_reader($name_sorted_sam_file);
		
		my $read_counter = 1;
		my $curr_read_name = undef;
		my %transcripts;
		
		
		$em_obj = EM->new($transcript_seqs_href);
		
		while (my $sam_entry = $sam_reader->get_next()) {
			
			my $read_name = $sam_entry->get_read_name();	
			my $trans_id = $sam_entry->get_scaffold_name();
			
			unless ($transcripts_passed_href->{$trans_id}) {
				next;
			}
			
			$transcripts{$trans_id} = 1;
			
			if ($curr_read_name && $curr_read_name ne $read_name) {
				$read_counter++;
				print STDERR "\r[$read_counter] reads read. ";
				if (%transcripts) {
					my @trans = keys %transcripts;
					# print "-adding read: @trans\n";
					$em_obj->add_read(@trans);
				}
				%transcripts = ();
			}
			
			$curr_read_name = $read_name;
		}
		
		
		

    
		# get last ones
		if (%transcripts) {
			$em_obj->add_read(keys %transcripts);
		}
		
		print STDERR "-done parsing reads, now running EM\n";
		
		$em_obj->run(max_iterations => $MAX_ITERATIONS_PHASE2,
					 min_delta_ML => $MIN_DELTA_ML);
		
		
		$em_obj->report_results();
		
		my @results = $em_obj->get_results();
		
		## check again to see if we need to filter some more low quality or poorly expressed transcripts
		$filtered_transcripts_flag = &examine_components_for_filtering(\@results, $transcripts_passed_href);
		

	}
	
	
	##########################
	####  Print Final Results
	##########################
	
	my $new_nameSorted_sam_file = "$OUT_PREFIX.nameSorted.sam";
    open (my $sam_ofh, ">$new_nameSorted_sam_file") or die "Error, cannot write to $new_nameSorted_sam_file";

	## write a new SAM file containing only alignments to the passed transcripts.
	my $sam_reader = new SAM_reader($name_sorted_sam_file);
	while (my $sam_entry = $sam_reader->get_next()) {
	
		my $read_name = $sam_entry->get_read_name();	
		my $trans_id = $sam_entry->get_scaffold_name();
	
		unless ($transcripts_passed_href->{$trans_id}) {
			next;
		}
	
        print $sam_ofh $sam_entry->toString() . "\n";

	}
	
	
	open (my $expr_ofh, ">$OUT_PREFIX.expression") or die "Error, cannot write to file $OUT_PREFIX.expression";
	open (my $fa_ofh, ">$OUT_PREFIX.fa") or die "Error, cannot write to file $OUT_PREFIX.fa";
	
	## write result file, order by component.
	my @results = $em_obj->get_results();
	@results = sort {$a->{trans_id} cmp $b->{trans_id}} @results;

	my $prev_component = undef;
	
	my $header = join("\t", "#transcript", "trans_length", "unique_map", "multi_map", "EM_frag_count", "FPKM", "%ComponentExpressed") . "\n";
	print $expr_ofh $header;
	print $header;
	
	## group by isoforms.
	my %component_to_isoforms;
	
	foreach my $result (@results) {
		
		
		
		my $trans_id = $result->{trans_id};
		$trans_id =~ /^(comp\d+)_/ or die "Error, cannot parse component acc: $trans_id";
		my $comp_id = $1;
		
		push (@{$component_to_isoforms{$comp_id}}, $result);
	}


	foreach my $component_id (keys %component_to_isoforms) {
		
		my @results = @{$component_to_isoforms{$component_id}};
		
		my $sum_fpkm = 0;
		foreach my $result (@results) {
			
			$sum_fpkm += $result->{FPKM};
		}

		foreach my $result (@results) {
			my $percent_gene_expressed = sprintf("%.2f", $result->{FPKM}/$sum_fpkm * 100);
					
			my $out_line = join("\t", 
								$result->{trans_id},
								$result->{length},
								$result->{unique_map},
								$result->{multi_map},
								$result->{expected_map},
								$result->{FPKM},
								"$percent_gene_expressed\%") . "\n";
			print $out_line;
			print $expr_ofh $out_line;
			
		
			my $trans_id = $result->{trans_id};
			my $transcript_seq = $transcript_seqs_href->{$trans_id};
			$transcript_seq =~ s/(\S{60})/$1\n/g;
			chomp $transcript_seq;
			
			print $fa_ofh ">" . $result->{trans_id} . " "
				. "length:" . $result->{length} . " "
				. "unique_map:" . $result->{unique_map} . " "
				. "multi_map:" . $result->{multi_map} . " "
				. "expected_map:" . $result->{expected_map} . " "
				. "FPKM:" . $result->{FPKM} . " "
				. "%ComponentExpressed:" . $percent_gene_expressed
				. "\n"
				. $transcript_seq . "\n";
			
		}
		print $expr_ofh "\n"; # spacer between components.
		print "\n"; 
	}
	
	close $expr_ofh;
	close $fa_ofh;

    ## create new coordSorted.sam file.
    my $cmd = "sort -T . -S 2G -k3,3 -k4,4n $new_nameSorted_sam_file > $OUT_PREFIX.coordSorted.sam";
    &process_cmd($cmd);
    
    
	return;
}


####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }
    return;
}

####
sub examine_components_for_filtering {
	my ($results_aref, $transcripts_passed_href) = @_;

	my %component_to_transcripts;
	
	my $filtered_flag = 0;

	foreach my $result (@$results_aref) {
		my $trans_id = $result->{trans_id};
		$trans_id =~ /^(comp\d+)_/ or die "Error, cannot parse component info from trans: $trans_id";
		my $comp_id = $1;
		push (@{$component_to_transcripts{$comp_id}}, $result);
	}

	foreach my $component_id (keys %component_to_transcripts) {
		
		my @results = @{$component_to_transcripts{$component_id}};
		
		my $sum_fpkm = 0;
		foreach my $result (@results) {
			my $fpkm = $result->{FPKM};
			$sum_fpkm += $fpkm;
		}
		
		foreach my $result (@results) {
			my $fpkm = $result->{FPKM};
			my $frags = $result->{expected_map};
			my $trans_id = $result->{trans_id};
			
			my $percent_component_expressed = 0;
			if ($sum_fpkm > 0) {
				$percent_component_expressed = $fpkm/$sum_fpkm*100;
			}
			
			if ($percent_component_expressed < $RETAIN_MIN_PERCENT_COMPONENT_EXPRESSED
				||
				$frags < 1
				)
			{
				$filtered_flag = 1;
				delete $transcripts_passed_href->{$trans_id};
			}
		}
	}

	return($filtered_flag);
}

####
sub require_unique_read_mapping {
	my ($read_to_trans_href) = @_;
	
	my $round = 0;

	while (1) {
		
		my %transcript_to_read_abundance;
		my %transcript_to_unique_read_count;
		
		foreach my $read (keys %$read_to_trans_href) {
			
			my $trans_href = $read_to_trans_href->{$read};
			my @trans = keys %$trans_href;
			
			if (scalar(@trans) == 1) {
				# uniquely mapped read
				$transcript_to_unique_read_count{$trans[0]}++;
				$transcript_to_read_abundance{$trans[0]}++;
			}
			else {
				foreach my $transcript (@trans) {
					$transcript_to_read_abundance{$transcript}++;
				}
			}

		}

		## see if there exists a transcript that lacks a uniquely placed read.
		my @non_unique_read_transcripts;
		foreach my $transcript (keys %transcript_to_read_abundance) {
			if (! exists $transcript_to_unique_read_count{$transcript}) {
				push (@non_unique_read_transcripts, $transcript);
			}
		}

		if (! @non_unique_read_transcripts) {
			last;
		}

		## remove the transcript with the lowest read count.
		@non_unique_read_transcripts  = sort {$transcript_to_read_abundance{$a}<=>$transcript_to_read_abundance{$b}} @non_unique_read_transcripts;
		
		if ($AGGRESSIVE) {
			&remove_transcript_from_read_mappings($read_to_trans_href, \@non_unique_read_transcripts);
			
			last;
		}
		else {
			$round++;
			print STDERR "\r[round $round] for removing non-uniquely mapped transcripts from component.     ";
			my $lowest_support_transcript = shift @non_unique_read_transcripts;
			&remove_transcript_from_read_mappings($read_to_trans_href, [$lowest_support_transcript]);
		}
		
	}
	
	#print STDERR "\n";
	
	return($read_to_trans_href);
}

####
sub remove_transcript_from_read_mappings {
	my ($read_to_trans_href, $trans_to_delete_aref) = @_;
	
	print STDERR " -removing transcripts: " . join(",", @$trans_to_delete_aref);
	
	my %to_delete = map { + $_ => 1 } @$trans_to_delete_aref;
	
	foreach my $read (keys %$read_to_trans_href) {
		
		my $trans_href = $read_to_trans_href->{$read};
		
		my @deletion;;
		foreach my $trans (keys %$trans_href) {
			if ($to_delete{$trans}) {
				push (@deletion, $trans);
			}
		}
		
		if (@deletion) {
			foreach my $trans (@deletion) {
				delete $trans_href->{$trans};
			}
		}
		
		if ($AGGRESSIVE) {
			# don't keep around empty reads, such as those that only mapped multiple times to transcripts now.
			unless (%$trans_href) {
				delete $read_to_trans_href->{$read};
			}
		}
	}
		

	return;
}
