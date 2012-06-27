#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use Nuc_translator;

my $usage = "usage: $0 454.fa\n\n";

my $roche454_reads = $ARGV[0] or die $usage;

my $pair_length = 250;
my $read_length = 76;

main: {

    my $fasta_reader = new Fasta_reader($roche454_reads);
    print STDERR "-parsing incoming 454 transcripts...";
    my %read_seqs = $fasta_reader->retrieve_all_seqs_hash();
    print STDERR "done.\n";
    
    my $num_454_trans = scalar (keys %read_seqs);
    my $counter_454 = 0;
    
    foreach my $read_acc (keys %read_seqs) {
        $counter_454++;
        print STDERR "\r[" . sprintf("%.2f%%  = $counter_454/$num_454_trans]     ", $counter_454/$num_454_trans*100);
        
        my $seq = $read_seqs{$read_acc};


        for (my $i = 0; $i <= length($seq); $i++) {
            
            
            my $left_read_seq = "";
            my $right_read_seq = "";
            my $ill_acc = $read_acc . "_p$i";
            

            my $left_start = $i - $read_length;
            if ($left_start >= 0) {
                $left_read_seq = substr($seq, $left_start, $read_length);
            }
                        
            my $right_start = $i;
            if ($right_start + $read_length -1 <= length($seq)) {
                $right_read_seq = substr($seq, $right_start, $read_length);
            }
            

            if ($left_read_seq) {
                print ">$ill_acc/1\n"
                    . "$left_read_seq\n";    
            }
            if ($right_read_seq) {
                print ">$ill_acc/2\n"
                    . "$right_read_seq\n";
            }
        }
    }
    

    print STDERR "\nDone.\n";
    
    exit(0);
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
sub capture_kmer_cov_text {
    my ($kmer_cov_file) = @_;
        
    my %kmer_cov;
    
    my $acc = "";
    open (my $fh, $kmer_cov_file) or die "Error, cannt open file $kmer_cov_file";
    while (<$fh>) {
        chomp;
        if (/>(\S+)/) {
            $acc = $1;
        }
        else {
            $kmer_cov{$acc} .= " $_";
        }
    }
    close $fh;

    return(%kmer_cov);
}


####
sub avg {
    my (@vals) = @_;

    if (scalar(@vals) == 1) {
        return($vals[0]);
    }
    

    my $sum = 0;
    foreach my $val (@vals) {
        $sum += $val;
    }
    
    my $avg = $sum / scalar(@vals);


    return(int($avg+0.5));
}

