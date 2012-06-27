#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

my $data_dir = "$FindBin::Bin/data";

my $usage = "usage: $0 KeggBlast.dump\n\n";

my $kegg_dump_file = $ARGV[0] or die $usage;

my $MIN_PER_ID = 0.50;
my $MIN_QUERY_COV = 0.50;


my %K_to_KOs = &parse_K_mappings("$data_dir/K-num_to_pathways.txt");
my %KO_descriptions = &parse_Kegg_pathway_descriptions("$data_dir/kegg_pathway_descriptions.txt");

my %best_hits;

open (my $fh, $kegg_dump_file) or die "Error, cannot open file $kegg_dump_file";
while (<$fh>) {
	chomp;
	my $line = $_;

	if (/^\#/) { next; }
	
	my ($trans_id, $gene_id, $locus, $gene_name, $q_descr, $q_sourceAcc, $q_sourceDB, $queryCoverage, $per_id, $score) = split(/\t/);
	
	unless ($per_id >= $MIN_PER_ID && $queryCoverage >= $MIN_QUERY_COV) {
		next;
	}
	
	my @Ks;
	while ($q_descr =~ /; (K\S+)/g) {
		my $K_num = $1;
		push (@Ks, $K_num);
	}

	if (@Ks) {
		my $struct = { 
			score => $score,
			K => \@Ks,
			line => $line,
		};
		
		my $prev_struct = $best_hits{$gene_id};
		
		if ( (! $prev_struct) || $prev_struct->{score} < $score) {

			$best_hits{$gene_id} = $struct;
		}
	}
}

foreach my $gene (keys %best_hits) {
	my $struct = $best_hits{$gene};
	
	my $line = $struct->{line};
	my ($trans_id, $gene_id, $locus, $gene_name, @rest) = split(/\t/, $line);
	
	my @Ks = @{$struct->{K}};

	foreach my $K (@Ks) {
		
		my @KOs = &get_KOs_via_K($K);
		foreach my $KO (@KOs) {
			my $descr = $KO_descriptions{$KO};
			print join("\t", $trans_id, $gene, $locus, $gene_name, $K, $KO, $descr) . "\n";
			
		}
	}
}


exit(0);

####
sub parse_K_mappings {
	my ($file) = @_;
	
	my %mappings;

	open (my $fh, $file) or die $!;
	while (<$fh>) {
		chomp;
		my ($k, $ko_list) = split(/\t/);
		my @KOs = split(/,/, $ko_list);
		
		$mappings{$k} = \@KOs;
	}
	close $fh;

	return(%mappings);
}


####
sub parse_Kegg_pathway_descriptions {
	my ($file) = @_;

	my %descr;

	open (my $fh, $file) or die "Error, cannot open file $file";
	while (<$fh>) {
		chomp;
		my ($ko, $descr) = split(/\t/);
		
		$descr{$ko} = $descr;
	}
	close $fh;

	return(%descr);
}

####
sub get_KOs_via_K {
	my ($K) = @_;

	my $KOs_aref = $K_to_KOs{$K};
	if (ref $KOs_aref) {
		return(@$KOs_aref);
	}
	else {
		return();
	}
}

