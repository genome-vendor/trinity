#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;

my $usage = "usage: $0 sampleA.RSEM.isoform.results sampleB.RSEM.isoform.results ...\n\n";

unless (@ARGV) {
    die $usage;
}

my @rsem_files = @ARGV;

unless (scalar @rsem_files > 1) {
    die $usage;
}

main: {


    my %data;
    
    foreach my $file (@rsem_files) {
        
        open (my $fh, $file) or die "Error, cannot open file $file";
        while (<$fh>) {
            chomp;
            my ($acc, $count, $rest) = split(/\t/);
            $data{$acc}->{$file} = $count;
        }
        close $fh;
    }

    my @filenames = @rsem_files;
    foreach my $file (@filenames) {
        $file = basename($file);
    }

    
    print join("\t", "", @filenames) . "\n";
    foreach my $acc (keys %data) {
        
        print "$acc";

        foreach my $file (@rsem_files) {

            my $count = $data{$acc}->{$file};
            unless (defined $count) {
                $count = "NA";
            }

            print "\t$count";
            
        }
        
        print "\n";
        
    }
    
    
    exit(0);
}
    
        
