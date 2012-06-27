#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);

use FindBin;

use Cwd;

$ENV{LC_ALL} = 'C';

my $util_dir = "$FindBin::Bin/../util";


my $usage = <<_EOUSAGE_;


#########################################################################################
#
# Required:
#
#  --iworm             inchworm assembled contigs
#
#  --left              left fragment file
#  --right             right fragment file
#
#  --seqType           fq|fa
#
# Optional (if strand-specific RNA-Seq):
#  
#  --SS_lib_type       RF or FR
#
#  --work_dir         directory to perform data processing (default: workdir.\$pid
#
#
###########################################################################################

_EOUSAGE_

	;


my $inchworm_contigs;
my $left_file;
my $right_file;
my $seqType;
my $SS_lib_type;
my $work_dir;

&GetOptions( 'iworm=s' => \$inchworm_contigs,
			 'left=s' => \$left_file,
			 'right=s' => \$right_file,
             'seqType=s' => \$seqType,
			 'SS_lib_type=s' => \$SS_lib_type,
			 'work_dir=s' => \$work_dir,
			 );


unless ($inchworm_contigs && $left_file && $right_file && $seqType) {
	die $usage;
}


unless ($work_dir) {
	$work_dir = "jaccard_clip_workdir";
}

main: {
	
	my $curr_dir = cwd();
	
	# create full paths to inputs, if not set.
	foreach my $file ($inchworm_contigs, $left_file, $right_file) {
		
		unless ($file =~ /^\//) {
			$file = "$curr_dir/$file";
		}
	}

	my $outdir = $work_dir;
	unless (-d $outdir) {
		mkdir ($outdir) or die "Error, cannot mkdir $outdir";
	}
	
    my $left_seqs = ($seqType eq 'fq') ? "left.fq" : "left.fa";
    my $right_seqs = ($seqType eq 'fq') ? "right.fq" : "right.fa";

	&process_cmd("ln -s $inchworm_contigs $outdir/iworm.fa") unless (-e "$outdir/iworm.fa");
	&process_cmd("ln -s $left_file $outdir/$left_seqs") unless (-e "$outdir/$left_seqs");
	&process_cmd("ln -s $right_file $outdir/$right_seqs") unless (-e "$outdir/$right_seqs");

	chdir $outdir or die "Error, cannot cd to $outdir";
	
   
    ## run the bowtie alignment pipeline
    
    my $cmd = "$util_dir/alignReads.pl --seqType $seqType --left $left_seqs --right $right_seqs --aligner bowtie --target iworm.fa ";
    if ($SS_lib_type) {
        $cmd .= "--SS_lib_type $SS_lib_type ";
    }
    &process_cmd($cmd);
    
		
	my $final_bam_file = ($SS_lib_type) ? "bowtie_out/bowtie_out.coordSorted.sam.+.bam" : "bowtie_out/bowtie_out.coordSorted.bam";
	
    my $alignment_file = "bowtie_alignments.sam";
    &process_cmd("samtools view $final_bam_file > $alignment_file");
    
    ## run Jaccard computation:
    my $jaccard_wig_file = "$alignment_file.J100.wig";
	my $frag_coords_file = "bowtie_alignments.sam.frag_coords";
    &process_cmd("$util_dir/SAM_ordered_pair_jaccard.pl --sam $alignment_file -W 100 > $jaccard_wig_file") unless (-s $jaccard_wig_file && -s $frag_coords_file);
    
    # The above creates a .frag_coords file.  Use this to compute fragment-level coverage
    my $frag_coverage_file = "$frag_coords_file.wig";
    &process_cmd("$util_dir/fragment_coverage_writer.pl $frag_coords_file > $frag_coverage_file");
    
    
	## define the transcript clip points:
    my $clips_file = "$jaccard_wig_file.clips";
	&process_cmd("$util_dir/jaccard_wig_clipper.pl --jaccard_wig $jaccard_wig_file --coverage_wig $frag_coverage_file  > $clips_file");

	## clip the inchworm transcripts:
	&process_cmd("$util_dir/jaccard_fasta_clipper.pl $inchworm_contigs  $clips_file > $inchworm_contigs.clipped.fa");
	
	
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
