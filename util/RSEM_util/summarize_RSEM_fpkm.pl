#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;
use File::Basename;

my $usage = <<__EOUSAGE__;

#################################################################################################
#
#  --transcripts <string>        transcripts fasta file (eg. Trinity.fasta)
#
#  --RSEM <string>               RSEM isoform abundance output file
#
#  --fragment_length <int>       mean length for a fragment targeted towards RNA-Seq (eg. 300)
#
# Optional:
#
#  --group_by_component          if Trinity assemblies, then groups data according to component
#
#################################################################################################


__EOUSAGE__

    ;

my $help_flag;
my $transcripts_fasta;
my $RSEM_file;
my $fragment_length;
my $group_by_component_flag = 0;


&GetOptions ( 'h' => \$help_flag,
              
              'transcripts=s' => \$transcripts_fasta,
              'RSEM=s' => \$RSEM_file,
              'fragment_length=i' => \$fragment_length,
              'group_by_component' => \$group_by_component_flag,
              
              );

my $MIN_EFF_LENGTH = 10;



unless ($transcripts_fasta && $RSEM_file && $fragment_length) {
    die $usage;
}

if ($help_flag) {
    die $usage;
}


main: {

    
    my %trans_lengths = &parse_trans_lengths($transcripts_fasta);

    my %data;
    
    my $sum_count = 0;

    open (my $fh, $RSEM_file) or die "Error, cannot open file $RSEM_file";
    while (<$fh>) {
        chomp;
        my ($trans, $count, $fraction) = split(/\t/);
        
        if ($fraction) {
            $fraction = sprintf("%.2e", $fraction);
        }
        else {
            $fraction = "NA";
        }
        
        my $comp = $trans;
        
        if ($group_by_component_flag) {
            $trans =~ /^(comp\d+_c\d+)/;
            $comp = $1 or die "Error, cannot extract trinity component name from accession: $trans";
        }

        my $length = $trans_lengths{$trans} or die "Error, cannot find length for transcript: $trans";


        $sum_count += $count;

        my $struct = { trans => $trans,
                       count => $count,
                       length => $length,
                       fraction => $fraction,
                   };
                       


        push (@{$data{$comp}}, $struct);
    }
    close $fh;


    ## compute the FPKM values
    foreach my $struct_list_aref (values %data) {

        foreach my $struct (@$struct_list_aref) {
            
            my $count = $struct->{count};
            my $length = $struct->{length};
            
            my $eff_length = $length - $fragment_length + 1;

            if ($eff_length < 1) {
                $eff_length = $MIN_EFF_LENGTH;
            }

            $struct->{eff_length} = $eff_length;
            
            my $fpkm = $count / ($eff_length/1e3) / ($sum_count / 1e6);

            $struct->{fpkm} = sprintf("%.2f", $fpkm);
        }
    }

    print "#Total fragments mapped to transcriptome: $sum_count\n";
    print join("\t", "transcript", "length", "eff_length", "count", "fraction", "fpkm", "\%comp_fpkm") . "\n";
    
    ## Report results
    foreach my $struct_list_aref (values %data) {
        
        my $sum_fpkm = 0;
        
        foreach my $struct (@$struct_list_aref) {
            $sum_fpkm += $struct->{fpkm};
        }

        # output
        foreach my $struct (@$struct_list_aref) {
            my $trans = $struct->{trans};
            my $count = $struct->{count};
            my $fraction = $struct->{fraction};
            my $length = $struct->{length};
            my $eff_length = $struct->{eff_length};
            
            my $fpkm = $struct->{fpkm};
            my $percent_comp_fpkm = 0;
            if ($sum_fpkm > 0) {
                $percent_comp_fpkm = sprintf("%.2f", $fpkm/$sum_fpkm*100);
            }
            
            print join("\t", $trans, $length, $eff_length, $count, $fraction, $fpkm, $percent_comp_fpkm) . "\n";
        }
        if ($group_by_component_flag) {
            print "\n"; # add spacer between components
        }
    }
    
    
    exit(0);

}

####
sub parse_trans_lengths {
    my ($trans_file) = @_;
    
    my $lengths_file = basename($trans_file) . ".seqLens";

    my %lengths;

    if (-s $lengths_file) {
        open (my $fh, $lengths_file) or die "Erorr, cannot read file $lengths_file";
        while (<$fh>) {
            chomp;
            my ($acc, $len) = split(/\t/);
            $lengths{$acc} = $len;
        }
        close $fh;
    
        return(%lengths);
    }
    else {
        
        open (my $ofh, ">$lengths_file") or die "Error, cannot write to $lengths_file";
        my $fasta_reader = new Fasta_reader($trans_file);
        while (my $seq_obj = $fasta_reader->next()) {

            my $acc = $seq_obj->get_accession();
            my $sequence = $seq_obj->get_sequence();

            my $seq_len = length($sequence);

            print $ofh "$acc\t$seq_len\n";
            
            $lengths{$acc} = $seq_len;
        }
       
        return(%lengths);
    }
}


        
  
