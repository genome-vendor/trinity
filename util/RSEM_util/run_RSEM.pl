#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

use Getopt::Long qw(:config no_ignore_case bundling);


my $usage = <<__EOUSAGE__;

#########################################################################
#
# --transcripts <string>           transcript fasta file
# 
# --name_sorted_bam <string>       name-sorted bam file  
#                                   eg. (sort -k1,1 -k3,3 file.sam > file.nameSorted.sam, then covert to bam using samtools)
#
# Optional:
# 
# --paired                         indicate that reads are paired
#
# --fragment_length <int>          for unpaired data:  indicate mean fragment length
#
# --SS_lib_type <string>           strand-specific library type:  paired('RF' or 'FR'), single('F' or 'R').
#
# --no_group_by_component          Trinity-mode, using 'components' as 'genes'
#
# --thread_count                   number of threads to use (default = 4)
#
# --debug                  retain intermediate files
#  
#########################################################################

__EOUSAGE__

    ;


my $help_flag;
my $transcripts;
my $bam_file;
my $paired_flag;
my $DEBUG_flag = 0;
my $fragment_length = 0;
my $SS_lib_type;
my $no_group_by_component = 0;
my $thread_count = 4;

&GetOptions ( 'h' => \$help_flag,
              'transcripts=s' => \$transcripts,
              'name_sorted_bam=s' => \$bam_file,
              'paired' => \$paired_flag,
              'debug' => \$DEBUG_flag,
              'fragment_length=i' => \$fragment_length,
              'SS_lib_type=s' => \$SS_lib_type,
              'no_group_by_component' => \$no_group_by_component,
              'thread_count=i' => \$thread_count,
              );



unless ($transcripts && $bam_file) {
    die $usage;
}


$paired_flag = ($paired_flag) ? "--paired-end" : "";


unless (defined $fragment_length) {
    die "Error, fragment length is not defined.\n";
}

if ($fragment_length && $paired_flag) {
    die "Error, do not specify --fragment_length with paired data";
}

if ($SS_lib_type) {
    unless ($SS_lib_type =~ /^(RF|FR|R|F)$/) {
        die "Error, do not recognize SS_lib_type: [$SS_lib_type]\n";
    }
    if ($paired_flag && length($SS_lib_type) != 2 ) {
        die "Error, SS_lib_type [$SS_lib_type] is not compatible with paired reads";
    }
}

if ( $thread_count !~ /^\d+$/ ) {
    die "Error, --thread_count value must be an integer";
}

my $RSEM_dir = "$FindBin::Bin/../../trinity-plugins/rsem";


main: {

    my $cmd = "$RSEM_dir/rsem-prepare-reference --no-polyA --no-bowtie";
    unless ($no_group_by_component) {
        my $trans_to_gene_map_file = &write_gene_to_trans_map_file($transcripts);
                
        $cmd .= " --transcript-to-gene-map $trans_to_gene_map_file";
    }
    $cmd .= " $transcripts TRANS";
    &process_cmd($cmd);
    
    my $keep_intermediate_files_opt = ($DEBUG_flag) ? "--keep-intermediate-files" : "";

    if ($fragment_length) {
        $fragment_length = "--fragment-length-mean $fragment_length";
    }
    else {
        $fragment_length = "";
    }

    my $SS_opt = "";
    if ($SS_lib_type) {
        if ($SS_lib_type =~ /^F/) {
            $SS_opt = "--forward-prob 1.0";
        }
        else {
            $SS_opt = "--forward-prob 0";
        }
    }
    
    $cmd = "$RSEM_dir/rsem-calculate-expression --no-qualities "
        . "$paired_flag "
        . "-p $thread_count "
        . "$fragment_length "
        . "$keep_intermediate_files_opt "
        . "$SS_opt "
        . "--bam $bam_file "
        . "TRANS "
        . "RSEM";
    
    &process_cmd($cmd);



    exit(0);
}


####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret: $ret";
    }
    
    return;
}


####
sub write_gene_to_trans_map_file {
    my ($transcripts_fasta_file) = @_;
    
        
    open (my $fh, $transcripts_fasta_file) or die "Error, cannot open file $transcripts_fasta_file";
    
    my $mapping_file = "$transcripts_fasta_file.component_to_trans_map";
    open (my $ofh, ">$mapping_file") or die "Error, cannot write to file: $mapping_file";
    
    while (<$fh>) {
        if (/>(comp\S+)/) {
            my $acc = $1;
            $acc =~ /^(comp\d+_c\d+)_seq\d+/ or die "Error, cannot parse the trinity component ID from $acc";
            my $comp_id = $1;
            print $ofh "$comp_id\t$acc\n";
        }
    }
    close $fh;
    close $ofh;

    return($mapping_file);
}
