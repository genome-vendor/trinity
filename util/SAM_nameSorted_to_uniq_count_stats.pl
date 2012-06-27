#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "\n\nusage: name_sorted_paired_reads.sam\n\n";

my $sam_file = $ARGV[0] or die $usage;

main: {


    my $prev_read_name = "";
    my $prev_scaff_name = "";

    my @reads;

    my %counts;

    my $sam_reader = new SAM_reader($sam_file);
    while ($sam_reader->has_next()) {
        

        my $read = $sam_reader->get_next();
        
        my $scaff_name = $read->get_scaffold_name();
        my $core_read_name = $read->get_core_read_name();
        
        if ($scaff_name ne $prev_scaff_name ||  $core_read_name ne $prev_read_name) {

            if (@reads) {
                &process_pairs(\@reads, \%counts);
                @reads = ();
            }
        }
        
        push (@reads, $read);
        

        $prev_read_name = $core_read_name;
        $prev_scaff_name = $scaff_name;
        
        
    }
    
    &process_pairs(\@reads, \%counts);
    

    my $sum_reads = 0;
    foreach my $count (values %counts) {
        $sum_reads += $count;
    }

    print "\n#read_type\tcount\tpct\n";
    foreach my $count_type (reverse sort {$counts{$a}<=>$counts{$b}} keys %counts) {
        my $count = $counts{$count_type};
        print "$count_type\t$count\t" . sprintf("%.2f", $count/$sum_reads*100) . "\n";
    }
    print "\n";
    
    exit(0);
}

####
sub process_pairs {
    my ($reads_aref, $counts_href) = @_;
    
    my @reads = @$reads_aref;
    
    my @left_reads;
    my @right_reads;

    foreach my $read (@reads) {
        if ($read->is_first_in_pair()) {
            push (@left_reads, $read);
        }
        elsif ($read->is_second_in_pair()) {
            push (@right_reads, $read);
        }
    }

    my $got_left_read = (@left_reads) ? 1 : 0;
    my $got_right_read = (@right_reads) ? 1: 0;
    

    ## check to see if we have proper pairs:
    my $got_propper_pair = 0;
    
  pair_search:
    foreach my $left_read (@left_reads) {
        
        my $aligned_pos = $left_read->get_aligned_position();
        
        
        foreach my $right_read (@right_reads) {
            
            if ($right_read->get_mate_scaffold_position() == $aligned_pos) {
                $got_propper_pair = 1;
                
                last pair_search;
            }
        }
    }
    
    if ($got_propper_pair) {
        $counts_href->{propper_pairs} += 2;
    }
    elsif ($got_left_read && $got_right_read) {
        $counts_href->{impropper_pairs} += 2;
    }
    elsif ($got_left_read) {
        $counts_href->{left}++;
    }
    elsif ($got_right_read) {
        $counts_href->{right}++;
    }
    else {
        die "Wassup....  not sure what I've got."; # shouldn't ever get here.
    }
    
    return;
}


    

        
    
