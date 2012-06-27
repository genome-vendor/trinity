#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);

use FindBin;

my $help_flag;

my $usage = <<__EOUSAGE__;

############################################################################################################################
#
# 
#  --target <string>         transcript fasta file (headers must have accessions that look like so:   gene_id;transcript_id
#  --query <string>          Trinity.fa
#  
#  --min_per_id <int>        min percent identity (default: 95)
#  --max_per_gap <int>  max percent gaps  (default: 5)
#
#  --forward_orient          strand-specific assemblies, only consider the forward strand.
#
#  --include_tiers           write the top tiers file
#
#  --allow_non_unique_mappings    a single query transcript counts towards all transcripts encapsulated  
#
#############################################################################################################################


__EOUSAGE__

    ;


my $target;
my $query;

my $min_per_id = 95;
my $max_per_gap = 5;

my $forward_orient = 0;
my $include_tiers = 0;

my $allow_non_unique_mappings = 0;

&GetOptions ( 'h' => \$help_flag,
              'target=s' => \$target,
              'query=s' => \$query,
              'min_per_id=i' => \$min_per_id,
              'forward_orient' => \$forward_orient,
              'max_per_gap=i' => \$max_per_gap,
              'include_tiers' => \$include_tiers,
              'allow_non_unique_mappings' => \$allow_non_unique_mappings,
              
    );


unless ($target && $query) {
    die $usage;
}

if ($help_flag) {
    die $usage;
}

my $util_dir = "$FindBin::Bin/util";


main: {

    my $cmd = "$util_dir/blat_full_length_mappings.pl --target $target --query $query "
        . " --min_per_id $min_per_id --max_per_gap $max_per_gap ";
    if ($forward_orient) {
        $cmd .= " --forward_orient ";
    }
    if ($allow_non_unique_mappings) {
        $cmd .= " --allow_non_unique_mappings ";
    }
    
    $cmd .= " > $query.pslx.maps";


    &process_cmd($cmd);  ## generates:  $query.pslx,  $query.pslx.FL_entries, $query.pslx.maps
    
    if ($include_tiers) {
        ## generate top tiers
        $cmd = "$util_dir/blat_top_tier_genes.pl $query.pslx > $query.pslx.tiers";
        &process_cmd($cmd);
    }
    
    ## generate summary statistics:
    $cmd = "$util_dir/blat_map_filter_with_isoforms.pl $query.pslx.maps | tee $query.pslx.maps.summary";
    &process_cmd($cmd);
    
    
    exit(0);

}

####
sub process_cmd {
    my ($cmd) = @_;

    print "CMD: $cmd\n";

    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}


