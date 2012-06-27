#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);

use FindBin;

use Cwd;
use File::Basename;

$ENV{LC_ALL} = 'C';

my $util_dir = "$FindBin::Bin/../util";

{
	my @progs_req = qw (bowtie-build bowtie);
	
	foreach my $prog (@progs_req) {
		unless (`which $prog` =~ /^\//) {
			die "Error, cannot find path to required program: $prog";
		}
	}
}

my $usage = <<_EOUSAGE_;

Note:  Uses bowtie to map reads.

#########################################################################################
#
# Required:
#
#  --target            fasta file containing target sequences to align to.
#
#  --left_fq           left fragment fq file
#  --right_fq          right fragment fq file
#
# Optional (if strand-specific RNA-Seq):
#  
#  --SS_lib_type       RF or FR
#
#  --work_dir         directory to perform data processing (default: bowtie.\$pid
#
#  -D                 max distance between pairs (default: 2000)
#  -C                 max number of read alignments reported (default: 1000)
#
###########################################################################################

_EOUSAGE_

	;


my $target_fasta;
my $left_fq;
my $right_fq;
my $SS_lib_type;
my $work_dir;

my $max_dist_between_pairs = 2000;
my $MAX_NUM_READ_ALIGNMENTS = 1000;

&GetOptions( 'target=s' => \$target_fasta,
			 'left_fq=s' => \$left_fq,
			 'right_fq=s' => \$right_fq,
			 'SS_lib_type=s' => \$SS_lib_type,
			 'work_dir=s' => \$work_dir,
             
             'D=i' => \$max_dist_between_pairs,
             'C=i' => \$MAX_NUM_READ_ALIGNMENTS,
             

			 );


unless ($target_fasta && $left_fq && $right_fq) {
	die $usage;
}


unless ($work_dir) {
	$work_dir = "bowtie.$$";
}

main: {
	
	my $curr_dir = cwd();
	
	# create full paths to inputs, if not set.
	foreach my $file ($target_fasta, $left_fq, $right_fq) {
		
		unless ($file =~ /^\//) {
			$file = "$curr_dir/$file";
		}
	}
	
	my $outdir = $work_dir;
	unless (-d $outdir) {
		mkdir ($outdir) or die "Error, cannot mkdir $outdir";
	}
	
	&process_cmd("ln -s $target_fasta $outdir/target.fa") unless (-e "$outdir/target.fa");
	&process_cmd("ln -s $left_fq $outdir/left.fq") unless (-e "$outdir/left.fq");
	&process_cmd("ln -s $right_fq $outdir/right.fq") unless (-e "$outdir/right.fq");

	chdir $outdir or die "Error, cannot cd to $outdir";
	
	&process_cmd("bowtie-build target.fa target") unless (-e "target.1.ebwt");
	
	## align left reads:
	&process_cmd("bowtie -q -S --sam-nohead -v 2 -k $MAX_NUM_READ_ALIGNMENTS target left.fq > left.fq.pre.sam ") unless (-s "left.fq.pre.sam");
    &process_cmd("$util_dir/SAM_filter_out_unmapped_reads.pl left.fq.pre.sam > left.fq.sam") unless (-s "left.fq.sam");
    &process_cmd("sort -T . -S 2G -k 1,1 -k 3,3 left.fq.sam > left.fq.nameSorted.sam") unless (-s "left.fq.nameSorted.sam");
    
	## align right reads:
    &process_cmd("bowtie -q -S -sam-nohead -v 2 -k $MAX_NUM_READ_ALIGNMENTS target right.fq >right.fq.pre.sam ") unless (-s "right.fq.pre.sam");
    &process_cmd("$util_dir/SAM_filter_out_unmapped_reads.pl right.fq.pre.sam > right.fq.sam") unless (-s "right.fq.sam");
	&process_cmd("sort -T . -S 2G -k 1,1 -k 3,3 right.fq.sam > right.fq.nameSorted.sam") unless (-s "right.fq.nameSorted.sam");
	
	## combine left and right into single file.
	&process_cmd("$util_dir/merge_left_right_nameSorted_SAMs.pl --left_sam left.fq.nameSorted.sam --right_sam right.fq.nameSorted.sam -D $max_dist_between_pairs -C $MAX_NUM_READ_ALIGNMENTS > combined.nameSorted.sam") unless (-s "combined.nameSorted.sam");
    
	## sort by coordinate.
	&process_cmd("sort -T . -S 2G -k 3,3 -k 4,4n combined.nameSorted.sam > combined.coordSorted.sam") unless (-s "combined.coordSorted.sam");
	
	
	my $final_coordSorted_sam_file = cwd() . "/combined.coordSorted.sam";
	my $final_nameSorted_sam_file = cwd() . "/combined.nameSorted.sam";

	if ($SS_lib_type) {
		## separate by strand:
		&process_cmd("$util_dir/SAM_strand_separator.pl combined.coordSorted.sam $SS_lib_type");
	
		$final_coordSorted_sam_file = cwd() . "/combined.coordSorted.sam.+.sam";
	
        &process_cmd("$util_dir/SAM_strand_separator.pl combined.nameSorted.sam $SS_lib_type");
        
        $final_nameSorted_sam_file = cwd() . "/combined.nameSorted.sam.+.sam";
        

    }
	
	
    my $final_coordSorted_sam_file_path = $curr_dir . "/" . basename($target_fasta) . ".bowtie.coordSorted.sam";
	my $final_nameSorted_sam_file_path = $curr_dir . "/" . basename($target_fasta) . ".bowtie.nameSorted.sam";
    
    my $final_coordSorted_bam_file_path = $curr_dir . "/" . basename($target_fasta) . ".bowtie.coordSorted.bam";
	    

	&process_cmd("ln -sf $final_coordSorted_sam_file $final_coordSorted_sam_file_path");
	&process_cmd("ln -sf $final_nameSorted_sam_file $final_nameSorted_sam_file_path");
    
	chdir ($curr_dir) or die "Error, cannot cd to $curr_dir";
	
    


	# convert to bam file.
	&process_cmd("samtools faidx $target_fasta");
	&process_cmd("samtools view -bt $target_fasta.fai $final_coordSorted_sam_file_path > $final_coordSorted_bam_file_path");
	
	
	
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

	return($ret);
}
