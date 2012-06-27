#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use FindBin;



my $usage = <<_EOUSAGE_;


##############################################################################################################
#
# Required:
#
#  --sam <string>                   sam-formatted file
#
# Optional:
#
#  --max_insert_size <int>          maximum insert size (default: 500) : used for pair coverage and jaccard
#  --min_insert_size <int>          minimum insert size (default: 100) : used for jaccard
#
#  --jaccard_win_length <int>       default: 100 (requires --jaccard)
#
##############################################################################################################



_EOUSAGE_

    ;


my $help_flag;

# required
my $sam_file;

# optional
my $max_insert_size = 500;
my $min_insert_size = 100;

my $jaccard_win_length = 100;


&GetOptions ( 'h' => \$help_flag,

              # required
              'sam=s' => \$sam_file,
              
              # optional
              'max_insert_size=i' => \$max_insert_size,
              'min_insert_size=i' => \$min_insert_size,
              
              'jaccard_win_length=i' => \$jaccard_win_length,
              
              );


if ($help_flag) {
    die $usage;
}


unless ($sam_file) {
    die $usage;
}


my $util_dir = "$FindBin::Bin";

main: {

        
    ## generate the paired coverage info:
    my $cmd = "$util_dir/SAM_to_frag_coords.pl --sam $sam_file "
        . "--max_insert_size $max_insert_size "
        . "--min_insert_size $min_insert_size ";  ## writes file: $sam_file.frag_coords
    
        
    &process_cmd($cmd) unless (-s "$sam_file.frag_coords");
    
    ## compute jaccard coeff info for pair-support
    $cmd = "$util_dir/ordered_fragment_coords_to_jaccard.pl --lend_sorted_frags $sam_file.frag_coords -W $jaccard_win_length ";

    &process_cmd($cmd);
    

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
