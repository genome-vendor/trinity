#!/usr/bin/env perl

use strict;
use warnings;
use File::Path;
use Carp;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");

use Nuc_translator;

my $usage = "usage: $0 partitions.gff alignments.sam [SS_lib_type=F,R,FR,RF]\n\n";

my $partitions = $ARGV[0] or die $usage;
my $alignments_sam = $ARGV[1] or die $usage;
my $SS_lib_type = $ARGV[2];

my $PARTS_PER_DIR = 100;


main: {


	my $partitions_dir = "Dir_${partitions}";
	mkdir ($partitions_dir) or die "Error, cannot mkdir $partitions_dir";

	open (my $track_fh, ">$partitions_dir.listing") or die $!;
	
	my %scaff_to_partitions = &parse_partitions($partitions);

	my @ordered_partitions;
	my $current_partition = undef;
	
	my $ofh;

	my $sam_ofh;
	
	my $partition_counter = 0;

	open (my $fh, $alignments_sam) or die "Error, cannot open file $alignments_sam";
	while (<$fh>) {
		#print;
		my $sam_line = $_;
		
		chomp;
		my @x = split(/\t/);
		my $acc = $x[0];
		my $scaff = $x[2];
		my $seq = $x[9];
		my $start = $x[3];
				
		unless (defined ($start) && $start =~ /^\d+$/ && defined($seq) && $seq =~ /^[gatcn]+$/i) {
			print "Line: $_ doesn't parse properly.  Skipping.\n";
			next;
		}

		if (! exists $scaff_to_partitions{$scaff}) { 
			# no partitions... should explore why this is.
			print STDERR "-warning, no read partitions defined for scaffold: $scaff\n";
			next; 
		} 
		
		my $flag = $x[1];
		
		my $aligned_strand = ($flag & 0x0010) ? '-' : '+';
		my $opposite_strand = ($aligned_strand eq '+') ? '-' : '+';
		
		if ($aligned_strand eq '-') {
			# restore to actual sequenced bases
			$seq = &reverse_complement($seq);
		}
		
		if ($SS_lib_type) {
			## got SS data
			
			my $transcribed_orient;

			my $pair_flag = ($flag & 0x0001) ? "PAIRED" : "notpaired";
			
			if ($pair_flag eq "notpaired") {
				if ($SS_lib_type !~ /^(F|R)$/) {
					confess "Error, read is not paired but SS_lib_type set to paired: $SS_lib_type\nread:\n$_";
				}
				
				if ($SS_lib_type eq "R") {
					$seq = &reverse_complement($seq);
				}
			}
			
			else {
				## Paired reads.
				if ($SS_lib_type !~ /^(FR|RF)$/) {
					confess "Error, read is paired but SS_lib_type set to unpaired: $SS_lib_type\nread:\n$_";
				}
						
				my $first_in_pair = $flag & 0x0040;
				if ( ($first_in_pair && $SS_lib_type eq "RF")
					 ||
					 ( (! $first_in_pair) && $SS_lib_type eq "FR")
					) {
					$seq = &reverse_complement($seq);
				}
			}
		}
		
	
		
		## prime ordered partitions if first entry or if switching scaffolds.
		if (! defined $current_partition || $scaff ne $current_partition->{scaff}) {
			@ordered_partitions = @{$scaff_to_partitions{$scaff}};
			$current_partition = shift @ordered_partitions;
			$partition_counter = 0; # reset
			close $ofh if $ofh;
			$ofh = undef;
		
			close $sam_ofh if $sam_ofh;
			$sam_ofh = undef;
			
		}
		
		## check to see if must advance partition
		if ($start > $current_partition->{rend}) {
			$current_partition = shift @ordered_partitions;
			$partition_counter++;
			close $ofh if $ofh;
			$ofh = undef;

			close $sam_ofh if $sam_ofh;
			$sam_ofh = undef;

			## note, might have now depleted partitions for this scaffold.
			

		}
		
		
		## check to see if we're in a partition
		if (defined($current_partition) && $start >= $current_partition->{lend} && $start <= $current_partition->{rend}) {
			# may need to start a new ofh for this partition if not already established.
			unless ($ofh) {
				# create new one.
				my $file_part_count = int($partition_counter/$PARTS_PER_DIR);
				my $outdir = "$partitions_dir/" . $current_partition->{scaff} . "/$file_part_count";
				
				mkpath($outdir) if (! -d $outdir);
				unless (-d $outdir) {
					die "Error, cannot mkdpath $outdir";
				}
				
				my $part_file = "$outdir/" . join("_", $current_partition->{lend}, $current_partition->{rend}) . ".reads";
				open ($ofh, ">$part_file") or die "Error, cannot write ot $part_file";
				print STDERR "-writing to $part_file\n";
				
				print $track_fh join("\t", $scaff, $current_partition->{lend}, $current_partition->{rend}, $part_file) . "\n";
				
				my $sam_part_file = "$outdir/" . join("_", $current_partition->{lend}, $current_partition->{rend}) . ".sam";
				open ($sam_ofh, ">$sam_part_file") or die "Error, cannot open $sam_part_file";
				

			}
			# write to partition
			print $ofh ">$acc\n$seq\n";
			print $sam_ofh $sam_line;
		}

	}
	close $fh;
	close $track_fh;
	
	close $ofh if $ofh;
	close $sam_ofh if $sam_ofh;

	exit(0);
}



####
sub parse_partitions {
	my ($partitions_file) = @_;
	
	my %scaff_to_parts;

	print STDERR "// parsing paritions.\n";
	my $counter = 0;

	open (my $fh, $partitions_file) or die "Error, cannot open file $partitions_file";
	while (<$fh>) {
		chomp;
		if (/^\#/) { next; }
		unless (/\w/) { next; }

		$counter++;
		print STDERR "\r[$counter]  ";
		
		my @x = split(/\t/);

		my $scaff = $x[0];
		my $lend = $x[3];
		my $rend = $x[4];
		my $orient = $x[6];
		
		push (@{$scaff_to_parts{$scaff}}, { scaff => $scaff,
											lend => $lend,
											rend => $rend, } );
		
	}
	
	close $fh;
	
	# should be sorted, but let's just be sure:
	foreach my $scaff (keys %scaff_to_parts) {
		@{$scaff_to_parts{$scaff}} = sort {$a->{lend}<=>$b->{lend}} @{$scaff_to_parts{$scaff}};
	}
	
	return(%scaff_to_parts);
}
			
	
