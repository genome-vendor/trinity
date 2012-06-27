#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");

use SAM_reader;
use SAM_entry;

use Getopt::Long qw(:config no_ignore_case bundling);

use Carp;
use Data::Dumper;

$ENV{LC_ALL} = 'C';  # critical for proper sorting using [system "sort -k1,1 ..."] within the perl script


my $usage = <<_EOUSAGE_;

#######################################################################
#
# Required:
#
# --sam <string>                SAM-formatted file
#
# Optional:
#
# --no_single                           no single reads
#
# --min_insert_size <int>               minimum distance between proper pair's starting positions. (default: 100)
# --max_insert_size <int>               maximum distance between proper pair's starting positions. (default: 500)
#
# -d                             debug mode
#
#######################################################################


_EOUSAGE_

    ;


my $sam_file;
my $full_flag = 0;
my $MAX_INSERT_SIZE = 500;
my $MIN_INSERT_SIZE = 100;
my $DEBUG = 0;
my $help_flag = 0;
my $NO_SINGLE;

&GetOptions( 
             'h' => \$help_flag,

             'sam=s' => \$sam_file,

             'no_single' => \$NO_SINGLE,
             'max_insert_size=i' => \$MAX_INSERT_SIZE,
             'min_insert_size=i' => \$MIN_INSERT_SIZE,
             
             'd' => \$DEBUG,

             );


if ($help_flag) {
    die $usage;
}


unless ($sam_file) {
    die $usage;
}


main: {
	

	my @frags;

    my $read_coords_file = "$sam_file.read_coords";

    &extract_read_coords($read_coords_file) unless (-s $read_coords_file);

    my $pair_coords_file = "$sam_file.frag_coords";
    &extract_frag_coords($read_coords_file, $pair_coords_file);
    
    
        
    
    exit(0);



}

####
sub extract_read_coords {
    my ($read_coords_file) = @_;

    ## Every read processed.

    print STDERR "-extracting read coordinates from $sam_file into $read_coords_file\n\n";
    
    my $sam_reader = new SAM_reader($sam_file);
    
    open (my $ofh, ">$read_coords_file") or die "Error, cannot write to file $read_coords_file";
    
	while ($sam_reader->has_next()) {
        
		my $sam_entry = $sam_reader->get_next();
		
        
		unless ($sam_entry->get_mate_scaffold_name() eq "=" || $sam_entry->get_mate_scaffold_name() eq $sam_entry->get_scaffold_name()) { next; }

    		
		my $scaffold = $sam_entry->get_scaffold_name();
        my $read_name = $sam_entry->get_core_read_name();
        my $full_read_name = $sam_entry->reconstruct_full_read_name();
        my $pair_side = ".";
        if ($full_read_name =~ m|/([12])$|) {
            $pair_side = $1;
        }
        

        my ($read_start, $read_end) = $sam_entry->get_genome_span();
        my $mate_scaff_pos = $sam_entry->get_mate_scaffold_position();
        
        
        print $ofh join("\t", $scaffold, $read_name, $pair_side, $read_start, $read_end) . "\n";
        
    }
    
    
    close $ofh;

    return;
        
}


####
sub extract_frag_coords {
    my ($read_coords_file, $pair_frag_coords_file) = @_;

    ## sort by scaffold, then by read name
    my $cmd = "sort -k1,1 -k2,2 -k4,4n $read_coords_file > $read_coords_file.sort_by_readname";
    &process_cmd($cmd);

    rename("$read_coords_file.sort_by_readname", $read_coords_file);
    
    ## define fragment pair coordinate span
    open (my $fh, $read_coords_file) or die $!;
    
    open (my $ofh, ">$pair_frag_coords_file") or die $!;
    
    
    my $first = <$fh>;
    chomp $first;
    while (my $second = <$fh>) {
        chomp $second;
        
        my ($scaffA, $readA, $readA_pair_side, $lendA, $rendA) = split(/\t/, $first);
        my ($scaffB, $readB, $readB_pair_side, $lendB, $rendB) = split(/\t/, $second);
        

        my $got_pair_flag = 0;
        if ($readA eq $readB && $scaffA eq $scaffB) {
            
            my @coords = sort {$a<=>$b} ($lendA, $rendA, $lendB, $rendB);
            my $min = shift @coords;
            my $max = pop @coords;
            
            my $insert_size = $max - $min + 1;
            if ($insert_size >= $MIN_INSERT_SIZE && $insert_size <= $MAX_INSERT_SIZE) {
                # treat as proper pair
                print $ofh join("\t", $scaffA, $readA, $min, $max) . "\n";
            }
            else {
                # treat as unpaired reads
                unless ($NO_SINGLE) {
                    print $ofh join("\t", $scaffA, $readA . "/$readA_pair_side", $lendA, $rendA) . "\n";
                    print $ofh join("\t", $scaffB, $readB . "/$readB_pair_side", $lendB, $rendB) . "\n";
                }
            }
    
            $first = <$fh>; # prime first
            chomp $first if $first;
        }
        else {
            # not paired
            unless ($NO_SINGLE) {
                print $ofh join("\t", $scaffA, $readA . "/$readA_pair_side", $lendA, $rendA) . "\n";
            }
            
            $first = $second;
            next;
        }
        
        

    }

    close $ofh;
    close $fh;


    $cmd = "sort -k1,1 -k3,3n $pair_frag_coords_file > $pair_frag_coords_file.coord_sorted";
    &process_cmd($cmd);

    rename("$pair_frag_coords_file.coord_sorted", $pair_frag_coords_file);
        
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

