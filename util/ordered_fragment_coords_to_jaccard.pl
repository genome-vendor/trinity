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
# --lend_sorted_frags <string>    lend-coordinate sorted fragment file
#
# Optional:
#
# -W <int>                       window length (default: 100)
# --pseudocounts <int>           default(1)           
# --full                         fully descriptive output format, rather than wig (default)
# --full_extreme                 includes fragment names and coordinates
# -v                             verbose, for progress-monitoring
#
# -d                             debug mode
#
#######################################################################




_EOUSAGE_

    ;


my $lend_sorted_frags_file;
my $window_length = 100;
my $full_flag = 0;
my $full_extreme_flag = 0;
my $VERBOSE = 0;
my $DEBUG = 0;
my $help_flag = 0;
my $pseudocounts = 1;

&GetOptions( 
             'h' => \$help_flag,
             
             'lend_sorted_frags=s' => \$lend_sorted_frags_file,
             "W=i" => \$window_length,
             
             'full' => \$full_flag,
             'full_extreme' => \$full_extreme_flag,
             
             'v' => \$VERBOSE,
             'd' => \$DEBUG,
             );

if ($help_flag) {
    die $usage;
}


unless ($lend_sorted_frags_file) {
    die $usage;
}

if ($full_extreme_flag) {
    $full_flag = 1;
}

main: {
	
        
    ## sort by rend
    my $rend_sorted_frags_file = "$lend_sorted_frags_file.sort_by_rend";
    my $cmd = "sort -k1,1 -k4,4n $lend_sorted_frags_file > $rend_sorted_frags_file";
    &process_cmd($cmd) unless (-s "$rend_sorted_frags_file");
    
    print STDERR "-processing jaccard pair sensor\n";
    &compute_jaccard_wig($lend_sorted_frags_file, $rend_sorted_frags_file, $window_length);
    
    
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
sub compute_jaccard_wig {
    my ($frags_sort_lend_file, $frags_sort_rend_file, $window_length) = @_;

    ## maintain two points, window length apart.
    ## at each point:  track accessions entering and leaving at each data point.


    ## Left window scanner
    my $left_window_scan_LEFT = Read_reader->new($frags_sort_lend_file);
    my $left_window_scan_RIGHT = Read_reader->new($frags_sort_rend_file);

    my $right_window_scan_LEFT = Read_reader->new($frags_sort_lend_file);
    my $right_window_scan_RIGHT = Read_reader->new($frags_sort_rend_file);
        
    my $curr_molecule = $left_window_scan_LEFT->get_curr_scaffold();
    print "variableStep chrom=$curr_molecule\n";
    
    my $window_lend = $left_window_scan_LEFT->get_curr_read_lend();
    my $window_rend = $window_lend + $window_length;
    
    my %frag_counter;

    my $num_single = 0;
    my $num_both = 0;
    
    my %extreme_verbose_frag_capture;

    while ($left_window_scan_LEFT->get_curr_line()) {
        
        ## see if we need to move on to the next molecule
        if ($left_window_scan_LEFT->get_curr_scaffold() gt $curr_molecule
            &&
            $left_window_scan_LEFT->get_curr_scaffold() eq $right_window_scan_LEFT->get_curr_scaffold()) {
            
            
            # init for next molecule
            $curr_molecule = $left_window_scan_LEFT->get_curr_scaffold();
            print "variableStep chrom=$curr_molecule\n";
            
            ## advance the other readers to ensure they're all synched up to the same molecule.
            $left_window_scan_RIGHT->advance_to_scaffold($curr_molecule);
            $right_window_scan_LEFT->advance_to_scaffold($curr_molecule);
            $right_window_scan_RIGHT->advance_to_scaffold($curr_molecule);
            

            %frag_counter = ();
            $num_single = 0;
            $num_both = 0;
            %extreme_verbose_frag_capture = ();
            

            $window_rend = $left_window_scan_LEFT->get_curr_read_lend();
            $window_lend = $window_rend - $window_length + 1;
            #if ($window_lend < 1) {
            #    $window_lend = 1;
            #}
        }
    

        {

            ## Fragment entering left
            ## process left window, lend of frag 
            
            my @acc_infos = $left_window_scan_LEFT->advance_get_accessions($curr_molecule, $window_lend, "LEND");
            foreach my $acc_info (@acc_infos) {
                my ($acc, $lend, $rend) = @$acc_info;
                
                

                if ($rend >= $window_lend) {
                    my $count = ++$frag_counter{$acc};
                    if ($count == 2) {
                        $num_single--;
                        $num_both++;
                    }
                    elsif ($count == 1) {
                        $num_single++;
                    }
                    if ($DEBUG) {
                        print "WINLEFT\tLEND\t$acc\t$window_lend\t$count\n";
                    }

                    $extreme_verbose_frag_capture{$acc} = "$acc\[$lend-$rend]" if ($full_extreme_flag);
                    
                
                }
            }
        }

        {

            ## Fragment exiting left
            ## process left window, rend of frag
            
            my @acc_infos = $left_window_scan_RIGHT->advance_get_accessions($curr_molecule, $window_lend-1, "REND");
            foreach my $acc_info (@acc_infos) {
                my ($acc, $lend, $rend) = @$acc_info;
                if ($rend == $window_lend-1) {
                    my $count = --$frag_counter{$acc};
                    if ($count == 1) {
                        $num_single++;
                        $num_both--;
                    }
                    elsif ($count == 0) {
                        $num_single--;
                        delete $frag_counter{$acc}; 
                    }
                    if ($DEBUG) {
                        print "WINLEFT\tREND\t$acc\t$window_lend\t$count\n";
                    }
                    $extreme_verbose_frag_capture{$acc} = "$acc\[$lend-$rend]" if ($full_extreme_flag);
                }
            }
        }

        {


            ## Fragment entering Right window
            ## process RIGHT window, lend of frag
            
            my @acc_infos = $right_window_scan_LEFT->advance_get_accessions($curr_molecule, $window_rend, "LEND");
            foreach my $acc_info (@acc_infos) {
                my ($acc, $lend, $rend) = @$acc_info;
                if ($rend >= $window_rend) {
                    my $count = ++$frag_counter{$acc};
                    if ($count == 2) {
                        $num_single--;
                        $num_both++;
                    }
                    elsif ($count == 1) {
                        $num_single++;
                    }
                    if ($DEBUG) {
                        print "WINRIGHT\tLEND\t$acc\t$window_rend\t$count\n";
                    }
                    $extreme_verbose_frag_capture{$acc} = "$acc\[$lend-$rend]" if ($full_extreme_flag);
                    
                }
            }
        }

        {

            ## Fragment exiting right window
            ## process right window, rend of frag
            
            my @acc_infos = $right_window_scan_RIGHT->advance_get_accessions($curr_molecule, $window_rend-1, "REND");
            foreach my $acc_info (@acc_infos) {
                my ($acc, $lend, $rend) = @$acc_info;
                if ($rend == $window_rend-1) {
                    my $count = --$frag_counter{$acc};
                    if ($count == 1) {
                        $num_single++;
                        $num_both--;
                    }
                    elsif ($count == 0) {
                        $num_single--;
                        delete $frag_counter{$acc}; 
                    }

                    if ($DEBUG) {
                        print "WINRIGHT\tREND\t$acc\t$window_rend\t$count\n";
                    }
                    $extreme_verbose_frag_capture{$acc} = "$acc\[$lend-$rend]" if ($full_extreme_flag);
                }
            }
        }
        
        
        ## compute jaccard coeff
        if ($window_lend >= 1) { # && ($num_single || $num_both)) {
            
            #my $jaccard = $num_both / ($num_single + $num_both);
            my $jaccard = ($num_both + $pseudocounts) / ($num_single + $num_both + $pseudocounts);
            
            $jaccard = sprintf("%.4f", $jaccard);
            my $mid = int( ($window_lend + $window_rend) / 2);
            
            if ($full_flag) {
                # include counts of single and paired frags
                print join("\t", $curr_molecule, $mid, $jaccard, "S:$num_single", "P:$num_both");
                if ($full_extreme_flag) {
                    my @paired;
                    my @single;
                    foreach my $frag (keys %frag_counter) {
                        if ($frag_counter{$frag} == 2) {
                            push (@paired, $frag);
                        }
                        elsif ($frag_counter{$frag} == 1) {
                            push (@single, $frag);
                        }
                        else {
                            die "Error, count not 1 or 2: " . Dumper(\%frag_counter);
                        }
                    }
                    
                    if ($num_single != scalar(@single) || $num_both != scalar(@paired)) {
                        die "Error, inconsistent counts: single: $num_single vs. " . scalar(@single) . ", and paired: $num_both vs. " . scalar(@paired);
                    }
                    
                    &print_frag_info(\@single, \@paired, \%extreme_verbose_frag_capture, $window_lend, $window_rend);
                    
                    
                }
                print "\n";
            }
            else {
                print join("\t", $mid, $jaccard) . "\n";
            }
        }
        
        
        ## slide window

        $window_rend++;
        $window_lend = $window_rend - $window_length;
        #if ($window_lend < 1) {
        #    $window_lend = 1;
        #}
        
        if ($VERBOSE && $window_lend % 1000 == 0) {
            print STDERR "\r[$window_lend] win lend  ";
        }
        
    }
    

    return;
    
}


####
sub print_frag_info {
    my ($single_aref, $paired_aref, $frag_info_href, $window_lend, $window_rend) = @_;
    
    my @left_single;
    my @right_single;
    
    foreach my $frag (@$single_aref) {
        my $frag_info = $frag_info_href->{$frag} or die "Error, no info for frag: $frag";
        
        $frag_info =~ /\[(\d+)-(\d+)\]$/ or die "Error, cannot parse coordinate info from $frag_info";
        my $frag_lend = $1;
        my $frag_rend = $2;
        
        if (&contains_pt($frag_lend, $frag_rend, $window_lend) && &contains_pt($frag_lend, $frag_rend, $window_rend)) {
            die "Error, frag $frag_info spans both edges of window[$window_lend-$window_rend] but was classified as single. ";
        }
        elsif (&contains_pt($frag_lend, $frag_rend, $window_lend) ) {
            push (@left_single, $frag_info);
        }
        elsif (&contains_pt($frag_lend, $frag_rend, $window_rend) ) {
            push (@right_single, $frag_info);
        }
        else {
            die "Error, single frag: $frag_info fails to overlap either edge of the window($window_lend-$window_rend)";
        }
        
    }
    
    my @paired_frags;
    foreach my $frag (@$paired_aref) {
        my $frag_info = $frag_info_href->{$frag} or die "Error, no info for frag: $frag";
        
        $frag_info =~ /\[(\d+)-(\d+)\]$/ or die "Error, cannot parse coordinate info from $frag_info";
        my $frag_lend = $1;
        my $frag_rend = $2;
        
        unless (&contains_pt($frag_lend, $frag_rend, $window_lend) && &contains_pt($frag_lend, $frag_rend, $window_rend)) {
            die "Error, paired frag $frag_info does NOT span both edges of window[$window_lend-$window_rend] but was classified as paired. ";
        }
        
        push (@paired_frags, $frag_info);
    }




    print "\tWINDOW:$window_lend-$window_rend"
        . "\tLEFT_SINGLE: " . join(",", @left_single)
        . "\tRIGHT_SINGLE: " . join(",", @right_single)
        . "\tPAIRS: " . join(",", @paired_frags);
    
    
    return;
    
}


####
sub contains_pt {
    my ($lend, $rend, $pt) = @_;

    if ($lend <= $pt && $pt <= $rend) {
        return(1);
    }
    else {
        return(0);
    }
}



##############################################################################
##############################################################################

package Read_reader;

use strict;
use warnings;
use Carp;

sub new {
    my $packagename = shift;
    my ($file) = @_;

    my $self = { file => $file,
                 fh => undef,
                 line => undef,
             };

    bless ($self, $packagename);
    
    $self->_init();

    return($self);
}
    

####
sub _init {
    my $self = shift;

    open (my $fh, $self->{file}) or die "Error, cannot open file: $self->{file}";
    $self->{fh} = $fh;

    $self->{line} = <$fh>;

    return;
}

####
sub advance_get_accessions {
    my $self = shift;
    my ($scaffold, $position, $type) = @_;
    

    #print STDERR "advance_get_accessions($scaffold, $position, $type)\n";

    unless ($type =~ /^LEND|REND$/) {
        confess "don't understand type: $type";
    }
    
    my @accs;
    while (1) {
        unless ($self->get_curr_line()) {
            return(@accs);
        }
        if ($scaffold ne $self->get_curr_scaffold()) {
            return(@accs);
        }

        my $acc = $self->get_curr_read_acc();
        my $lend = $self->get_curr_read_lend();
        my $rend = $self->get_curr_read_rend();
        
        my $pos = ($type eq "LEND") ? $lend : $rend;

        
        #print "$type\t$pos\n";
        
        if ($pos <= $position) {
            
            unless ($acc =~ m|/[12]$|) { # ignore individual reads, only consider complete fragments
                push (@accs, [$acc, $lend, $rend]);
            }
            
            $self->advance_line();
        }
        else {
            return(@accs);
        }
        

    }

    croak "shouldn't ever get here.";
    
}


####
sub get_curr_line {
    my $self = shift;
    
    return($self->{line});
}

####
sub get_curr_fields {
    my $self = shift;
    
    my $line = $self->get_curr_line();
    chomp $line;
    my @x = split(/\t/, $line);
    return(@x);
}

####
sub get_curr_scaffold {
    my $self = shift;
    my @x = $self->get_curr_fields();
    
    return($x[0]);
}

sub advance_to_scaffold {
    my $self = shift;
    my ($scaff) = @_;

    my $curr_scaff = $self->get_curr_scaffold();
    while (defined($curr_scaff) && $curr_scaff ne $scaff) {
        $self->advance_line();
        $curr_scaff = $self->get_curr_scaffold();    
    }
    
    return;
}


####
sub get_curr_read_acc {
    my $self = shift;
    
    my @x = $self->get_curr_fields();
    
    return($x[1]);
}


####
sub get_curr_read_lend {
    my $self = shift;
    my @x = $self->get_curr_fields();
    return($x[2]);
}

####
sub get_curr_read_rend {
    my $self = shift;
    my @x = $self->get_curr_fields();
    return($x[3]);
}

####
sub advance_line {
    my $self = shift;
    my $fh = $self->{fh};
    my $next_line = <$fh>;
    
    while ($next_line && $next_line =~ m|/[12]$|) {
        $next_line = <$fh>;
    }
        
    
    $self->{line} = $next_line;
    
    return($next_line);
}
    

