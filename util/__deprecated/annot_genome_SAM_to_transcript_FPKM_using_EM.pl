#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib/");
use Gene_obj;
use GFF3_utils;
use GTF_utils;
use SAM_reader;
use SAM_entry;
use EM;
use Overlap_info;
use Getopt::Long qw(:config no_ignore_case bundling);
use Data::Dumper;

my $usage = <<_EOUSAGE_;

#########################################################################
#
#  ### Annot settings:
#
#  --gff3 <string>     gff3 file name
# 
#    OR
#
#  --gtf <string>      gtf file name 
#
#    AND
#
# ### SAM settings:
#
#  --coord_sorted_sam <string>  genome coordinate-sorted sam file name
#
# ### rna-seq info
#
#  --frag_length      length of an RNA-Seq fragment (** not the length of a read, but rather the mean length of a fragment from which the read was derived **)
#
# # Below are optional:
#  
#  --SS_lib_type       one of [RF,FR,F,R,XS].  if set to XS, then using the XS:A:[+-] attribute for assignment.
#   
#  --FUZZY             require compatible overlap of read/transcript but do not require containment of entire read.
#  --REQUIRE_PAIRINGS  transcript must have both paired fragment ends mapped in order to be counted.
#
#  --max_iterations    for EM algorithm, default 1000
#  --min_delta_ML      for EM algorithm, default 0.01
#
#  --outfile           output file name (default:  sam_filename.fpkm)
#
###################################################################################################


_EOUSAGE_

	;




my $genes_gff3;
my $genes_gtf;
my $sam_file;
my $help_flag;
my $SS_lib_type;
my $VERBOSE = 0;
my $min_delta_ML = 0.01;
my $max_iterations = 1000;
my $FUZZY_OVERLAP = 0;
my $REQUIRE_PAIRINGS = 0;
my $outfile = "";
my $frag_length;;

&GetOptions ( 'h' => \$help_flag,
			  'gff3=s' => \$genes_gff3,
			  'gtf=s' => \$genes_gtf,
			  'coord_sorted_sam=s' => \$sam_file,
              'SS_lib_type=s' => \$SS_lib_type,
              'max_iterations=i' => \$max_iterations,
              'min_delta_ML=f' => \$min_delta_ML,
			  'v' => \$VERBOSE,
              'FUZZY' => \$FUZZY_OVERLAP,
              'REQUIRE_PAIRINGS' => \$REQUIRE_PAIRINGS,
              'outfile=s' => \$outfile,
              'frag_length=i' => \$frag_length,
              );


if ($help_flag) { die $usage; }

unless ($sam_file) {
    die $usage;
}

unless ($genes_gff3 || $genes_gtf) { 
	die $usage;
}


unless ($frag_length) {
    die $usage;
}

unless ($outfile) {
    $outfile = "$sam_file.fpkm";
}


main: {
	
    my $gene_obj_indexer_href = {};
	my $contig_to_gene_list_href;

    print STDERR "-parsing gene annotations file\n";
	if ($genes_gff3) {
	
		## associate gene identifiers with contig id's.
		$contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($genes_gff3, $gene_obj_indexer_href);
	}
	else {
		# GTF mode
		$contig_to_gene_list_href = &GTF_utils::index_GTF_gene_objs_from_GTF($genes_gtf, $gene_obj_indexer_href);
	}
    
    
    print STDERR "-organizing transcript features.\n";
    
    ## Populate Genomic Features
	
	my %chr_to_features;
	
    my %acc_to_seqs;

    my $total_transcripts = 0;

	foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
		my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
		foreach my $gene_id (@gene_ids) {
			my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
			
			unless (ref $gene_obj_ref) {
				die "Error, no gene_obj for gene_id: $gene_id";
			}
			
			my $strand = $gene_obj_ref->get_orientation();

			
			my $scaffold = $asmbl_id;
            
			my $max_isoform_cdna_length = 0;
			
			foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
                
                $total_transcripts++;
                
                my $isoform_id = join("::", $gene_id, $isoform->{Model_feat_name}); ## embed the gene identifier in with the transcript id, so we can tease it apart later.
				my $isoform_length = 0;
				
                my @coordset = &get_coordset_for_isoform($isoform);
                my ($lend, $rend) = ($coordset[0]->[0], $coordset[$#coordset]->[1]);
                
                my $length = &sum_coordset_segments(\@coordset);
                $acc_to_seqs{$isoform_id} = 'N' x $length;
                
                push (@{$chr_to_features{$scaffold}}, { 
                    acc => $isoform_id,
                    name => $isoform->{com_name},
                    lend => $lend,
                    rend => $rend,
                    coordset => [@coordset],
                    strand => $strand,
                    length => $length,
                      } );
                
            }
		}
	}
    
    
    ## Check to see if we just want to generate a summary
    if (-s $outfile) {
        print STDERR "Output file: $outfile already exists.  Should we simply use this for generating the summary? ";
        my $response = <STDIN>;
        if ($response =~ /^y/i) {
            &generate_expressed_transcript_report($outfile, $total_transcripts);
            exit(0);
        }
    }
        
    
	##
	## Map reads to genes
	##

	my %gene_id_to_reads_mapped;
	
	my %fragment_tracker;
	my $read_counter = 0;
	my %acc_to_read_count;

	
	## for current processing.
	my $curr_scaff = "";
	my @chr_features;	
	my @container;  ## holds current features.

    
	
    print STDERR "-examining read alignments, mapping to annotated features.\n";
    

    my $trans_mapping_file = "read_trans_mappings.$$.txt";
    open (my $ofh, ">$trans_mapping_file") or die $!;
    
    
    my $sam_reader = new SAM_reader($sam_file);
    while ($sam_reader->has_next()) {
        
        my $sam_entry = $sam_reader->get_next();
                
        if ($sam_entry->is_query_unmapped()) { 
            next;
        }
        
        $read_counter++;
        
        if ($read_counter % 1000 == 0) {
            print STDERR "\r[$read_counter] reads examined.     ";
        }
        
        my $read_acc =  $sam_entry->get_core_read_name();
        my $full_read_acc = $sam_entry->reconstruct_full_read_name();
                
        my $molecule = $sam_entry->get_scaffold_name();
        my $position_lend = $sam_entry->get_aligned_position();
        
        my ($read_align_coords_aref, $read_coords_aref) = $sam_entry->get_alignment_coords();
        
        my $orient = "?";
        if ($SS_lib_type) {
            if ($SS_lib_type eq "XS") {
                ## examine XS:A attribute
                my $sam_line = $sam_entry->toString();
                if ($sam_line =~ /XS:A:([\+\-])/) {
                    $orient = $1;
                }
            }
            else {
                $orient = $sam_entry->get_query_transcribed_strand($SS_lib_type);
            }
        }
        
        
        $fragment_tracker{$read_acc} = 1;
      			
        
        if ($molecule ne $curr_scaff) {
            ## reset
            if (exists $chr_to_features{$molecule}) {
                
                @chr_features = sort {$a->{lend}<=>$b->{lend}} @{$chr_to_features{$molecule}};
            }
            else {
                print STDERR "Warning, no chr features for molecule: $molecule\n";
            }
            @container = (); # clear
            $curr_scaff = $molecule;
        }
        
        @container = sort {$a->{rend}<=>$b->{rend}} @container;
        ## purge current contained features that no longer overlap position.
        while (@container && $container[0]->{rend} < $position_lend) {
            shift @container;
        }
        
        ## collect features that overlap 
        while (@chr_features && 
               $chr_features[0]->{lend} <= $position_lend) {
            
            my $feature = shift @chr_features;
            if ($feature->{rend} >= $position_lend) {
                
                push (@container, $feature);
            }
            else {
                # no op, feature gets tossed.
            }
        }
        
        if (@container) {
            
            foreach my $feature (@container) {
                my $acc = $feature->{acc};
                if ( ($FUZZY_OVERLAP && &Overlap_info::compatible_overlap($feature->{coordset}, $read_align_coords_aref) )
                      ||
                      (&Overlap_info::compatible_overlap_A_contains_B($feature->{coordset}, $read_align_coords_aref))
                      ) {
                    
                    print STDERR "\r[$read_counter] reads examined.     Adding read count to $acc  " if $VERBOSE;
                    print $ofh join("\t", $full_read_acc, $acc) . "\n";
                    
                    
                }
                
            }
            
            
            
            
        }
        else {
            #print STDERR "\nNo annot mapped to $_\n";
            
        }
        
        
        #if ($read_counter > 100000) { last; } ###DEBUG
    }
    
    
    close $ofh;

    

    ## Do EM, estimate abundance
    my $num_fragments_mapped = scalar (keys %fragment_tracker);
    %fragment_tracker = (); # no longer needed.


    &process_cmd("sort -k1,1 $trans_mapping_file > $trans_mapping_file.sort");
    
    my $em_obj = new EM(\%acc_to_seqs, $frag_length);
    
    $em_obj->{_total_reads} = $num_fragments_mapped;  # use total mapped to genome rather than total mapped to transcripts.
    
    my %seen_trans;
    my %read_trans_mapping;
    
    my $prev_core_read = "";
    
    open (my $fh, "$trans_mapping_file.sort") or die $!;
    while (<$fh>) {
        chomp;
        my ($full_read_name, $trans) = split(/\t/);
        my $core_read_name = $full_read_name;
        $core_read_name =~ s|/\d$||;
                
        
        if ($prev_core_read ne $core_read_name) {
            ## process current read set
            my @candidate_trans = keys %seen_trans;
            my @mapped_trans = ();
            if ($REQUIRE_PAIRINGS) {
                foreach my $candidate (@candidate_trans) {
                    if ($read_trans_mapping{"${prev_core_read}/1"}->{$candidate}
                        &&
                        $read_trans_mapping{"${prev_core_read}/2"}->{$candidate}) {
                        ## got pairing
                        push (@mapped_trans, $candidate);
                    }
                }
            }
            else {
                @mapped_trans = @candidate_trans;
            }
            
            
            if (@mapped_trans) {
                $em_obj->add_read(@mapped_trans);
            }
            
            # reinit
            %seen_trans = ();
            %read_trans_mapping = ();
        }
        
        $seen_trans{$trans} = 1;
        $read_trans_mapping{$full_read_name}->{$trans} = 1;
        
        $prev_core_read = $core_read_name;
        
    }
    
    if (%seen_trans) {
        
        ## get last ones   TODO: make this a function, since this is just code duplication from above.
        
        my @candidate_trans = keys %seen_trans;
        my @mapped_trans = ();
        if ($REQUIRE_PAIRINGS) {
            foreach my $candidate (@candidate_trans) {
                if ($read_trans_mapping{"${prev_core_read}/1"}->{$candidate}
                    &&
                    $read_trans_mapping{"${prev_core_read}/2"}->{$candidate}) {
                    ## got pairing
                    push (@mapped_trans, $candidate);
                }
            }
        }
        else {
            @mapped_trans = @candidate_trans;
        }
        
        
        if (@mapped_trans) {
            $em_obj->add_read(@mapped_trans);
        }
        
        
    }


    ## Run EM to estimate FPKM values.
    $em_obj->run(max_iterations => $max_iterations, min_delta_ML => $min_delta_ML, verbose => $VERBOSE);

    
    ## Report output
    open ($ofh, ">$outfile") or die "Error, cannot write to $outfile";


    print $ofh "# Total fragments mapped to genome: $num_fragments_mapped\n";
    print $ofh $em_obj->report_results();
    close $ofh;
        
    unlink ($trans_mapping_file, "$trans_mapping_file.sort");
    
    &generate_expressed_transcript_report($outfile, $total_transcripts);
    

    exit(0);
}


####
sub generate_expressed_transcript_report {
    my ($fpkm_file, $total_transcripts) = @_;

    my $Rscript = "$fpkm_file.R";
    open (my $ofh, ">$Rscript");
    print $ofh "source(\"$FindBin::Bin/R/expression_analysis_lib.R\")\n";
    print $ofh "png(\"$fpkm_file.genes_vs_minFPKM.png\", width=1000, height=500)\n";
    print $ofh "plot_expressed_gene_counts(\"$fpkm_file\", total=$total_transcripts, title=\"expressed transcripts vs. min FPKM\", fpkm_range=seq(0,5,0.01), outfile=\"$fpkm_file.genes_vs_minFPKM.dat\")\n";
    print $ofh "dev.off()\n";
    close $ofh;
    
    system("R --vanilla -q < $Rscript");
    
    return;
}
    




####
sub get_coordset_for_isoform {
    my ($isoform) = @_;

    my @coordsets;
    
    my @exons = $isoform->get_exons();
						
    foreach my $exon (@exons) {
        my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();

        push (@coordsets, [$lend, $rend]);
    }

    @coordsets = sort {$a->[0]<=>$b->[0]} @coordsets;


    return(@coordsets);
}


####
sub sum_coordset_segments {
    my ($coordset_aref) = @_;

    my $sum = 0;
    foreach my $coordset (@$coordset_aref) {
        
        my ($lend, $rend) = @$coordset;
        
        $sum += $rend - $lend + 1;
    }

    return($sum);
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
