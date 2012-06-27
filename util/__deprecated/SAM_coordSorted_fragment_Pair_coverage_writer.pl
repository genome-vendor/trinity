#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "usage: $0 coordSorted.file.sam [max_dist=500]\n\n";

my $sam_file = $ARGV[0] or die $usage;
my $MAX_DIST = $ARGV[1] || 500; # maximum distance between alignment starts of paired reads.

main: {
	
	my %scaffold_to_coverage; # will retain all coverage information.
	
	my $sam_reader = new SAM_reader($sam_file);
	
	my $current_scaff = undef;
	
	my $counter = 0;
	while ($sam_reader->has_next()) {

		$counter++;
		if ($counter % 1000 == 0) {
			print STDERR "\r[$counter]  ";
		}

		my $sam_entry = $sam_reader->get_next();

        if ($sam_entry->is_query_unmapped()) { next; }
        
        unless ($sam_entry->is_proper_pair()) { next; }
        
		
		if (%scaffold_to_coverage && $sam_entry->get_scaffold_name() ne $current_scaff) {
			&report_coverage(\%scaffold_to_coverage);
			%scaffold_to_coverage = ();
		}
		
		$current_scaff = $sam_entry->get_scaffold_name();

		&add_coverage($sam_entry, \%scaffold_to_coverage);
		

	}

	if (%scaffold_to_coverage) {
		&report_coverage(\%scaffold_to_coverage);
	}
		

	exit(0);

}


####
sub report_coverage {
	my ($scaffold_to_coverage_href) = @_;
	
	## output the coverage information:
	foreach my $scaffold (sort keys %$scaffold_to_coverage_href) {
		
		print "variableStep chrom=$scaffold\n";
		
		my @coverage = @{$scaffold_to_coverage_href->{$scaffold}};
		
		for (my $i = 1; $i <= $#coverage; $i++) {
			my $cov = $coverage[$i] || 0;
			
			print "$i\t$cov\n";
		}
		
	}
	
	return;
}


####
sub add_coverage {
	my ($sam_entry, $scaffold_to_coverage_href) = @_;
	
	
	## If the entry is paired and the query alignment
	## is less than the paired position, then fill in the gap.
		
	my $scaffold = $sam_entry->get_scaffold_name();
	my $scaff_pos = $sam_entry->get_scaffold_position();
	
	my @genome_coords;
	
    my ($read_lend, $read_rend) = $sam_entry->get_genome_span();
    my $mate_pos = $sam_entry->get_mate_scaffold_position();

    if ($mate_pos > $read_lend) {
        $read_rend = $mate_pos - 1;
    }

    if ($read_rend - $read_lend + 1 > $MAX_DIST) {
        return;
    }
    
    ## add coverage:
    for (my $i = $read_lend; $i <= $read_rend; $i++) {
        $scaffold_to_coverage_href->{$scaffold}->[$i]++;
    }
		

	return;
}
	
