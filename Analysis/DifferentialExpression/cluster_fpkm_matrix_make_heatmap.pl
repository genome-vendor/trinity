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
#  --matrix <string>     matrix.RAW.normalized.FPKM
#
#  Optional:
#
#  --output <string>     prefix for output file (default: "\${matrix_file}.heatmap")
#
#  --log2_median_center      
#
####################################################################################


__EOUSAGE__

    ;


my $matrix_file;
my $output_prefix = "";
my $LOG2_MEDIAN_CENTER = 0;

&GetOptions (  'matrix=s' => \$matrix_file,
               'output=s' => \$output_prefix,
               'log2_median_center' => \$LOG2_MEDIAN_CENTER,
               );

unless ($matrix_file) {
    die $usage;
}

unless ($output_prefix) {
    $output_prefix = "$matrix_file.heatmap";
}


main: {
    
    
    my $R_script = "$output_prefix.R";
    
    open (my $ofh, ">$R_script") or die "Error, cannot write to $R_script";
    
    print $ofh "library(cluster)\n";
    print $ofh "library(gplots)\n";
    print $ofh "library(Biobase)\n";
    
    print $ofh "data = read.table(\"$matrix_file\", header=T, com=\'\', sep=\"\\t\")\n";
    print $ofh "rownames(data) = data[,1] # set rownames to gene identifiers\n";
    print $ofh "data = data[,2:length(data[1,])] # remove the gene column since its now the rowname value\n";
    print $ofh "data = as.matrix(data) # convert to matrix\n";
    if ($LOG2_MEDIAN_CENTER) {
        print $ofh "data = log2(data+1)\n";
        print $ofh "data = t(scale(t(data), scale=F)) # center rows, mean substracted\n";
            ;
    }
    #print $ofh "g_dist = dist(1-cor(t(data)), method='euclidean')\n";

    print $ofh "hc_genes = agnes(data, diss=FALSE, metric=\"euclidean\") # cluster genes\n";
    #print $ofh "hc_genes = hclust(as.dist(1-cor(t(data), method=\"spearman\")), method=\"complete\") # cluster conditions\n";
    #    ;
    print $ofh "hc_samples = hclust(as.dist(1-cor(data, method=\"spearman\")), method=\"complete\") # cluster conditions\n";
        ;
    print $ofh "myheatcol = redgreen(75)\n";
    print $ofh "gene_partition_assignments <- cutree(as.hclust(hc_genes), k=6);\n";
    print $ofh "partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)\n";
    print $ofh "gene_colors = partition_colors[gene_partition_assignments]\n";
    
    print $ofh "save(list=ls(all=TRUE), file=\"all.RData\")\n";
    
    print $ofh "postscript(file=\"$output_prefix.eps\", horizontal=FALSE, width=8, height=18, paper=\"special\");\n";
    print $ofh "heatmap.2(data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors, scale=\"none\", density.info=\"none\", trace=\"none\", key=TRUE, keysize=1.2, cexCol=2.5, margins=c(8,8), lhei=c(0.3,2), lwid=c(2.5,4))\n";
    print $ofh "dev.off()\n";
    
    
    close $ofh;
    
    &process_cmd("R --vanilla -q < $R_script");


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

