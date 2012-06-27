#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use Getopt::Long qw(:config no_ignore_case bundling);

my $usage = <<_EOUSAGE_;

########################################################################################################
#
#  Required:
#
#  --coord_sorted_SAM <string>      coordinate-sorted SAM file.
#
#  -J  <int>                       value corresponding to the 75th percentile of intron lengths (coverage partitions within this length are joined)
#  -I  <int>                       maximum intron length  (reads with longer intron lengths are ignored)
#                             
#  *If Strand-specific, specify:
#  --SS_lib_type <string>          library type:  if single: F or R,  if paired:  FR or RF
#
#
#  Optional:
#
#  --jaccard_clip                 mitigates false fusions of genes.  ** requires paired-reads. **
#
#
########################################################################################################


_EOUSAGE_

	;

my $help_flag;

my $partition_join_size;
my $max_intron_length;
my $SAM_file;
my $SS_lib_type = "";
my $jaccard_clip_flag = 0;

&GetOptions ( 'h' => \$help_flag,

			  'coord_sorted_SAM=s' => \$SAM_file,
			  'SS_lib_type=s' => \$SS_lib_type,
			  
			  'J=i' => \$partition_join_size,
			  'I=i' => \$max_intron_length,
			  			  
			  'jaccard_clip' => \$jaccard_clip_flag,
	);


if ($help_flag) {
	die $usage;
}

unless (
	$SAM_file
	&& $partition_join_size
	&& $max_intron_length
		) {
	die $usage;
}

if ($SS_lib_type && $SS_lib_type !~ /^(F|R|FR|RF)$/) {
	die "Error, invalid --SS_lib_type, only F, R, FR, or RF are possible values";
}

my $UTIL_DIR = "$FindBin::Bin/";

main: {

	my @sam_info;

	if ($SS_lib_type) {
		my ($plus_strand_sam, $minus_strand_sam) = ("$SAM_file.+.sam", "$SAM_file.-.sam");
		if (-s $plus_strand_sam && $minus_strand_sam) {
			print STDERR "-strand partitioned SAM files already exist, so using them instead of re-creating them.\n";
		}
		else {
			my $cmd = "$UTIL_DIR/SAM_strand_separator.pl $SAM_file $SS_lib_type";
			&process_cmd($cmd);
		}

		push (@sam_info, [$plus_strand_sam, '+'], [$minus_strand_sam, '-']);
	}
	else {
		push (@sam_info, [$SAM_file, '+']);
	}
			

	foreach my $sam_info_aref (@sam_info) {
				
		my ($sam, $strand) = @$sam_info_aref;
		
		## compute base coverage:
		my $cmd = "$UTIL_DIR/SAM_coordSorted_fragment_coverage_writer2.pl $sam $max_intron_length > $sam.wig";
		&process_cmd($cmd) unless (-s "$sam.wig");

		## define partitions based on coverage:
		$cmd = "$UTIL_DIR/define_SAM_coverage_partitions2.pl $sam.wig $strand $max_intron_length > $sam.partitions";
		&process_cmd($cmd) unless (-s "$sam.partitions");

		## join partitions within join size:
		$cmd = "$UTIL_DIR/join_partitions_within_range.pl $sam.partitions $partition_join_size > $sam.partitions.join$partition_join_size.gff";
		&process_cmd($cmd) unless (-s "$sam.partitions.join$partition_join_size.gff");

		my $partitions_gff = "$sam.partitions.join$partition_join_size.gff";

		if ($jaccard_clip_flag) {

			$cmd = "$UTIL_DIR/SAM_ordered_pair_jaccard.pl --sam $sam > $sam.jaccard";
			&process_cmd($cmd) unless (-s "$sam.jaccard");

			$cmd = "$UTIL_DIR/jaccard_wig_clipper.pl --jaccard_wig $sam.jaccard > $sam.jaccard.clips";
			&process_cmd($cmd) unless (-s "$sam.jaccard.clips");

			$cmd = "$UTIL_DIR/segment_GFF_partitions.pl $sam.partitions.join$partition_join_size.gff $sam.jaccard.clips > jac.$strand.partitions.gff";
			&process_cmd($cmd) unless (-s "jac.$strand.partitions.gff");

			$partitions_gff = "jac.$strand.partitions.gff";
		}
		


		## extract reads per partition
		$cmd = "$UTIL_DIR/extract_reads_per_partition.pl $partitions_gff $sam $SS_lib_type";
		&process_cmd($cmd) unless (-d "Dir_jac.$strand.partitions.gff");
		
		
		
	}

	print "##\nDone\n##\n\n";


	exit(0);
	
	

}


####
sub process_cmd {
	my ($cmd) = @_;

	print STDERR "CMD: $cmd\n";
	
	my $ret = system($cmd);

	if ($ret) {
		die "Error, command $cmd died with ret $ret";
	}

	return;
}
