#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use Cwd;
use File::Basename;

use Getopt::Long qw(:config no_ignore_case bundling);


## check for required programs
{

	my @required_progs = qw (blat psl2sam.pl samtools);
	foreach my $prog (@required_progs) {
		my $path = `which $prog`;
		unless ($path =~ /^\//) {
			die "Error, path to required $prog cannot be found";
		}
	}
}




$ENV{LC_ALL} = 'C';

my $usage = <<_EOUSAGE_;

################################################################################################################
#
#  --left and --right    (if paired reads)
#     or
#  --single              (if unpaired reads)
#
#  Required inputs:
#   
#  --genome            multi-fasta file containing the genome sequences (should be named {refName}.fa )
#
#  --seqType          fa | fq    (fastA or fastQ format)
#
# Optional:
#
#  --SS_lib_type      strand-specific library type:  single: F or R  paired: FR or RF
#                                examples:  single RNA-Ligation method:  F
#                                           single dUTP method: R
#                                           paired dUTP method: RF
#
#  -I    maximum intron length  (default: 10000);
#
#  -o    output directory
#
#  --trim_short_terminal_segments     (trim off short terminal alignment segments that are mostly noise. Default: 10) 
# 
#  -P   min percent identity based on full sequence length  (default: 95)
# 
#  --blat_top_hits  (default: 20 in paired mode, 1 in single mode)
#
#  -C  final top hits reported  (default: 1)  (only applies to paired mode)
#
#  If paired mode:
#
#     --max_dist_between_pairs             default (2000) 
#
####################################################################################################################


_EOUSAGE_

	;


my $help_flag;
my $genome_db; 
my $left_file;
my $right_file;
my $single_file;
my $max_intron = 10000;
my $output_directory = "output_$$";
my $trim_short_terminal_segment_length = 10;
my $min_per_ID = 95;
my $BLAT_top_hits = 20;
my $MIN_TOP_HITS = 1;
my $max_dist_between_pairs = 2000;
my $seqType;
my $SS_lib_type;

&GetOptions ( 'h' => \$help_flag,
			  
			  ## required inputs
			  'left=s' => \$left_file,
			  'right=s' => \$right_file,
			  
			  'single=s' => \$single_file,
			  

			  "genome=s" => \$genome_db,
			  "seqType=s" => \$seqType,
			  

			  ## Optional:
			  "SS_lib_type=s" => \$SS_lib_type,
			  
			  "I=i" => \$max_intron,
			  
			  'o=s' => \$output_directory,
			  
			  'trim_short_terminal_segments=i' => \$trim_short_terminal_segment_length,

			  'P=i' => \$min_per_ID,
			  
			  'blat_top_hits=i' => \$BLAT_top_hits,
			  
			  'C=i' => \$MIN_TOP_HITS,
			  
			  'max_dist_between_pairs=i' => \$max_dist_between_pairs,
			  
	);


if ($help_flag) { die $usage; }

unless ($genome_db && -s $genome_db) { 
	die $usage;
}

unless ($seqType && $seqType =~ /^(fq|fa)$/) {
	die $usage;
}

unless ( ($single_file && -s $single_file)
		 || 
		 ($left_file && -s $left_file
		  && $right_file && -s $right_file)) {
	die $usage;
}


if ($single_file) {
	$BLAT_top_hits = $MIN_TOP_HITS;
}

if ($SS_lib_type && $SS_lib_type !~ /^(F|R|FR|RF)$/) {
	die "Error, SS_lib_type must be one of the following: (F, R, FR, RF)  ";
}


main: {


	my $start_dir = cwd();

	$single_file = &build_full_paths($single_file, $start_dir) if $single_file;
	$left_file = &build_full_paths($left_file, $start_dir) if $left_file;
	$right_file = &build_full_paths($right_file, $start_dir) if $right_file;
	$genome_db = &build_full_paths($genome_db, $start_dir);
	
	my $work_dir;
	if ($output_directory =~ /^\//) {
		$work_dir = $output_directory;
	}
	else {
		$work_dir = "$start_dir/$output_directory";
	}
	
	&process_cmd("mkdir -p $work_dir") unless (-d $work_dir);
	
	my $util_dir = "$FindBin::Bin/../util";
	
	my @entries;
	
	if ($single_file) {
		push (@entries, ["single_dir", "single_fa", $single_file]);
	}
	else {
		push (@entries, 
			  ["left_dir", "left_fa", $left_file],
			  ["right_dir", "right_fa", $right_file]);
	}

	unless ($genome_db =~ /^\//) {
		## full path not given, must reconstruct it.
		$genome_db = "$work_dir/$genome_db";
	}
	
	unless (-s "genome.fa") {
		
		# prep the genome here for converting sam to bam later.
		
		my $cmd = "ln -s $genome_db genome.fa";
		&process_cmd($cmd);
	}
	
	unless (-s "genome.fa.fai") {
		my $cmd = "samtools faidx genome.fa";
		&process_cmd($cmd);
	}
	

	foreach my $fa_info (@entries) {
		
		## Prep work area:
		chdir $work_dir or die "Error, cannot cd to $work_dir";
		
		unless (-s "genome.fa") {
			
			# prep the genome here for converting sam to bam later.
			
			my $cmd = "ln -s $genome_db genome.fa";
			&process_cmd($cmd);
			
			$cmd = "samtools faidx genome.fa";
			&process_cmd($cmd);
		}


		my ($target_dir, $target_fa, $trans_file) = @$fa_info;
		
		unless (-d $target_dir) {
			mkdir $target_dir or die "Error, cannot mkdir $target_dir";
		}
		
		# prep target dir with symlinks to components needed.
		if ($seqType eq "fq") {
			my $cmd = "$util_dir/fastQ_to_fastA.pl -I $trans_file > $target_dir/$target_fa";
			&process_cmd($cmd) unless (-e "$target_dir/$target_fa");
		}
		else {
			## already fasta format:
			my $cmd = "ln -s $trans_file $target_dir/$target_fa";
			&process_cmd($cmd);
		}
		
		my $cmd = "ln -s $genome_db $target_dir/genome.fa";
		&process_cmd($cmd) unless (-e "$target_dir/genome.fa");
		

		## Work in Target_dir
		chdir ($target_dir) or die "Error, cannot cd to $target_dir";
		
		## Prep sequences>
		if ($seqType eq "fq") {
			$cmd = "$util_dir/fastQ_to_tab.pl -I $trans_file > $target_fa.tab";
			&process_cmd($cmd) unless (-s "$target_fa.tab");
		}
		else {
			$cmd = "$util_dir/fasta_to_tab.pl < $trans_file > $target_fa.tab";
			&process_cmd($cmd) unless (-s "$target_fa.tab");
		}
		
		$cmd = "sort -T . -S 2G -k 1,1 -k 3,3 $target_fa.tab > $target_fa.sorted.tab";
		&process_cmd($cmd) unless (-s "$target_fa.sorted.tab");

		## run blat
		
		$cmd = "$util_dir/run_BLAT_shortReads.pl genome.fa $target_fa $max_intron $target_fa.psl";
		&process_cmd($cmd) unless (-s "$target_fa.psl");
				
		## convert to sam
		$cmd = "psl2sam.pl -q 0 -r 0 $target_fa.psl > $target_fa.psl.sam";
		&process_cmd($cmd) unless (-s "$target_fa.psl.sam");
		
		## sort by name
		$cmd = "sort -T . -S 2G -k 1,1 -k 3,3 $target_fa.psl.sam > $target_fa.psl.nameSorted.sam";
		&process_cmd($cmd) unless (-s "$target_fa.psl.nameSorted.sam");
		
		## add sequences to sam file.
		$cmd = "$util_dir/blat_sam_add_reads2.pl $target_fa.psl.nameSorted.sam $target_fa.sorted.tab > $target_fa.psl.nameSorted.wReads.sam";
		&process_cmd($cmd) unless (-s "$target_fa.psl.nameSorted.wReads.sam");
		
		## capture top hits
		$cmd = "$util_dir/top_blat_sam_extractor.pl $target_fa.psl.nameSorted.wReads.sam $BLAT_top_hits $min_per_ID > $target_fa.nameSorted.sam";
		&process_cmd($cmd) unless (-s "$target_fa.nameSorted.sam");
		
		
	}
	

	chdir $work_dir or die "Error, cannot cd to $work_dir";
	
	## merge into single sam file, setting flags properly

	my $outfile_basename = basename($output_directory);
	
	if ($single_file) {
		
		my $cmd = "sort -T . -S 2G -k 3,3 -k 4,4n single_dir/single_fa.nameSorted.sam > $outfile_basename.coordSorted.sam";
		&process_cmd($cmd);
        
    		
	}
	else {
		## paired mode:
	
		
		my $cmd = "$util_dir/merge_left_right_nameSorted_SAMs.pl --left_sam left_dir/left_fa.nameSorted.sam --right_sam right_dir/right_fa.nameSorted.sam -C $MIN_TOP_HITS -D $max_dist_between_pairs > combined.nameSorted.sam";
		&process_cmd($cmd) unless (-s "combined.nameSorted.sam");
		
		## sort by coordinate.
		
		$cmd = "sort -T . -S 2G -k 3,3 -k 4,4n combined.nameSorted.sam > $outfile_basename.coordSorted.sam";
		&process_cmd($cmd);
	}

	
	# report splice junctions and remove short terminal exons that are more likely noise.
	my $cmd = "$FindBin::Bin/../Inchworm/bin/cigar_tweaker $outfile_basename.coordSorted.sam genome.fa $trim_short_terminal_segment_length | sort -T . -S 2G -k 3,3 -k 4,4n >  $outfile_basename.coordSorted.spliceAdjust.sam";
	&process_cmd($cmd) unless (-s "$outfile_basename.coordSorted.spliceAdjust.sam");
	
	
	# add transcribed orientation info:
	if ($SS_lib_type) {
		$cmd = "$util_dir/SAM_set_transcribed_orient_info.pl $outfile_basename.coordSorted.spliceAdjust.sam $SS_lib_type > $outfile_basename.coordSorted.spliceAdjust.senseOrientSet.sam";
		&process_cmd($cmd) unless (-s "$outfile_basename.coordSorted.spliceAdjust.senseOrientSet.sam");
	}
	else {
		# not strand-specific, keep as is and don't disrupt current flow (so use expected output name)
		$cmd = "ln -s  $work_dir/$outfile_basename.coordSorted.spliceAdjust.sam $outfile_basename.coordSorted.spliceAdjust.senseOrientSet.sam";
		&process_cmd($cmd);
	}
	
	# convert to bam format
	$cmd = "samtools view -bt genome.fa.fai -S $outfile_basename.coordSorted.spliceAdjust.senseOrientSet.sam > $outfile_basename.coordSorted.bam";
	&process_cmd($cmd) unless (-s "$outfile_basename.coordSorted.bam");
	
	rename("$outfile_basename.coordSorted.spliceAdjust.senseOrientSet.sam", "$start_dir/$outfile_basename.coordSorted.sam") or die "Error, cannot relocate $outfile_basename.coordSorted.spliceAdjust.senseOrientSet.sam to $start_dir";
	
	rename("$outfile_basename.coordSorted.bam", "$start_dir/$outfile_basename.coordSorted.bam") or die "Error, cannot relocate $outfile_basename.coordSorted.bam to $start_dir";
	
	chdir($start_dir) or die "Error, cannot cd back to $start_dir";

    ## provide name-sorted SAM
    $cmd = "sort -T . -S 2G -k 1,1 -k 3,3 $outfile_basename.coordSorted.sam > $outfile_basename.nameSorted.sam";
    &process_cmd($cmd);
	
	$cmd = "samtools index $outfile_basename.coordSorted.bam";
	&process_cmd($cmd);

	if ($SS_lib_type) {
		## strand-specific
		## separate the sam based on strand, and create separate bam files.  (for convenience sake)
		
		$cmd = "$util_dir/SAM_strand_separator.pl $outfile_basename.coordSorted.sam $SS_lib_type";
		&process_cmd($cmd);

        $cmd = "$util_dir/SAM_strand_separator.pl $outfile_basename.nameSorted.sam $SS_lib_type";
        &process_cmd($cmd);
        
		foreach my $sam_file ("$outfile_basename.coordSorted.sam.+.sam", "$outfile_basename.coordSorted.sam.-.sam") {
			
			if (-s $sam_file) {

				my $bam_file = $sam_file;
				$bam_file =~ s/\.sam$/\.bam/;
				
				$cmd = "samtools view -bt genome.fa.fai $sam_file > $bam_file";
				&process_cmd($cmd);
				
				$cmd = "samtools index $bam_file";
				&process_cmd($cmd);
			}
		}
	}
	
	
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
sub build_full_paths {
	my ($path, $start_dir) = @_;
	
	if ($path && $path !~ /^\//) {
		$path = "$start_dir/$path";
	}

	return($path);
}
