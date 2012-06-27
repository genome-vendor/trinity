#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);


my $usage = <<__EOUSAGE__;

#################################################################################### 
#
# Required:
#
#  --matrix    matrix.RAW.normalized.FPKM
#
# Optional:
#
#  -P          p-value cutoff for FDR  (default: 0.001)
#
#  -C          min abs(log2(a/b)) fold change (default: 2  (meaning 2^2 or 4-fold).
#
#  --output    prefix for output file (default: "diffExpr.P\${Pvalue}_C\${C})
#
#  --DESeq     parse \*DESeq.output files instead of edgeR results.txt files. 
#
####################################################################################


__EOUSAGE__

    ;


my $matrix_file;
my $p_value = 0.001;
my $log2_fold_change = 2;
my $output_prefix = "";
my $FORCE_FLAG = 0;
my $help_flag = 0;
my $DESeq_mode = 0;

&GetOptions (  'h' => \$help_flag,
               
               'matrix=s' => \$matrix_file,
               'P=f' => \$p_value,
               'C=f' => \$log2_fold_change,
               'output=s' => \$output_prefix,
               'FORCE' => \$FORCE_FLAG, # for exploratory purposes.
               'DESeq' => \$DESeq_mode,
               );


if ($help_flag) {
    die $usage;
}

unless ($matrix_file) {
    die $usage;
}

unless ($output_prefix) {
    $output_prefix = "diffExpr.P${p_value}_C${log2_fold_change}";
}


main: {

    ## get list of genes that meet the criterion:
    my %diffExpr;
    if ($DESeq_mode) {
        %diffExpr = &parse_DESeq_results($p_value, $log2_fold_change);
    }
    else {
        %diffExpr = &parse_result_files_find_diffExp($p_value, $log2_fold_change);
    }
    
    unless (%diffExpr) {
        die "Error, no differentially expressed transcripts identified at cuttoffs: P:$p_value, C:$log2_fold_change";
    }

    my $diff_expr_matrix = "$output_prefix.matrix";
    open (my $ofh, ">$diff_expr_matrix") or die "Error, cannot write to file $diff_expr_matrix";
    
    open (my $fh, $matrix_file) or die "Error, cannot read file $matrix_file";
    my $header = <$fh>;
    print $ofh $header;
    my $count_diff_expr = 0;
    while (<$fh>) {
        my $line = $_;
        my @x = split(/\t/);
        my $acc = $x[0];
        
        if ($diffExpr{$acc}) {
            print $ofh $line;
            $count_diff_expr++;
        }
    }

    close $ofh;

    unless ($FORCE_FLAG) {
        unless ($count_diff_expr == scalar(keys %diffExpr)) {
            die "Error, number of diff expr genes identified in $matrix_file is $count_diff_expr, and != to " . scalar(keys %diffExpr) . " number of genes identified by edgeR as diffExpressed.";
        }
    }
    
    &cluster_diff_expressed_transcripts($diff_expr_matrix);
    

    exit(0);
    
}

####
sub cluster_diff_expressed_transcripts {
    my ($diff_expr_matrix_file) = @_;

    my $R_script = "$diff_expr_matrix_file.R";
    
    open (my $ofh, ">$R_script") or die "Error, cannot write to $R_script";
    
    print $ofh "library(cluster)\n";
    print $ofh "library(gplots)\n";
    print $ofh "library(Biobase)\n";
    
    print $ofh "data = read.table(\"$diff_expr_matrix_file\", header=T, com=\'\', sep=\"\\t\")\n";
    print $ofh "rownames(data) = data[,1] # set rownames to gene identifiers\n";
    print $ofh "data = data[,2:length(data[1,])] # remove the gene column since its now the rowname value\n";
    print $ofh "data = as.matrix(data) # convert to matrix\n";
    print $ofh "data = log2(data+1)\n";
    print $ofh "centered_data = t(scale(t(data), scale=F)) # center rows, mean substracted\n";
    print $ofh "hc_genes = agnes(centered_data, diss=FALSE, metric=\"euclidean\") # cluster genes\n";
    print $ofh "hc_samples = hclust(as.dist(1-cor(centered_data, method=\"spearman\")), method=\"complete\") # cluster conditions\n";
    print $ofh "myheatcol = redgreen(75)\n";
    print $ofh "gene_partition_assignments <- cutree(as.hclust(hc_genes), k=6);\n";
    print $ofh "partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)\n";
    print $ofh "gene_colors = partition_colors[gene_partition_assignments]\n";
    
    print $ofh "save(list=ls(all=TRUE), file=\"${R_script}.all.RData\")\n";
    
    print $ofh "postscript(file=\"$diff_expr_matrix_file.heatmap.eps\", horizontal=FALSE, width=8, height=18, paper=\"special\");\n";
    print $ofh "heatmap.2(centered_data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors, scale=\"none\", density.info=\"none\", trace=\"none\", key=TRUE, keysize=1.2, cexCol=2.5, margins=c(15,15), lhei=c(0.3,2), lwid=c(2.5,4))\n";
    print $ofh "dev.off()\n";
    
    
    close $ofh;
    
    &process_cmd("R --vanilla -q < $R_script");

    return;
    

=notes from zehua regarding heatmap.2

you need to change the margin of the heatmap in the command heatmap.2:

margins=c(8,8),

The first number is the margin on the bottom (column name), and the second number is for the margin on the right (i.e., the row names). So you can increase the second number and you should get larger space for row names.

If you want to use a column order as you specified, then you can just turn off the ordering of the column by setting the following options in the heatmap.2 command:

dendrogram='row', Colv=FALSE

This will allow you to have a heatmap with column order from the data you provided.

=cut



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
sub parse_result_files_find_diffExp {
    my ($max_p_value, $min_abs_log2_fold_change) = @_;

    my %diff_expr;

    my @result_files = <*results.txt>;
    
    unless (@result_files) {
        die "Error, cannot find edgeR \*results.txt files. Be sure to run this from within the edgeR output directory";
    }
    

    foreach my $file (@result_files) {
        open (my $fh, $file) or die "Error, cannot open file $file";
        my $header = <$fh>;
        while (<$fh>) {
            my ($gene, $log_conc, $log_fold_change, $p_value, $fdr) = split(/\s+/);
            
            if (abs($log_fold_change) >= $min_abs_log2_fold_change
                &&
                $fdr <= $max_p_value) {

                $diff_expr{$gene} = 1;
            }
        }
        close $fh;
    }

    return(%diff_expr);
}

####
sub parse_DESeq_results {
    my ($p_value, $log2_fold_change) = @_;

    my %diff_expr;

    my @files = <*.DESeq.output>;

    unless (@files) {
        die "Error, no \*.DESeq.output files identified in current directory.";
    }
    
    foreach my $file (@files) {
        open (my $fh, $file) or die "Error, cannot open file $file";
        my $header = <$fh>;
        while (<$fh>) {
            my @x = split(/\t/);
            
            my $feature = $x[0];
            my $log2fc = $x[5];
            my $adjpval = $x[7];
            
            if ( $log2fc ne "NA" && abs($log2fc) > $log2_fold_change
                 &&
                 $adjpval ne "NA" && $adjpval <= $p_value) {

                $diff_expr{$feature} = 1;
            }
        }
        close $fh;
    }
 
    unless (%diff_expr) {
        die "Error, no diff expressed genes identified";
    }
                
    return(%diff_expr);
}
