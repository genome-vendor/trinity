#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;

my $usage = <<__EOUSAGE__;

###################################################################################
#
# -K <int>     define K clusters
# -R <string>  the filename for the store RData (file.all.RData)
#
###################################################################################


__EOUSAGE__

    ;


my $K;
my $help_flag = 0;
my $R_data_file;

&GetOptions ( 'h' => \$help_flag,
              'K=i' => \$K,
              'R=s' => \$R_data_file,
              );


if ($help_flag) {
    die $usage;
}

unless ($K && $R_data_file) {
    die $usage;
}



main: {
    
    unless (-s $R_data_file) {
        die "Error, cannot find pre-existing R-session data as file: $R_data_file";
    }
    
    
    my $R_script = "__tmp_define_${K}_clusters.R";
    
    open (my $ofh, ">$R_script") or die "Error, cannot write to file $R_script";

    print $ofh "library(cluster)\n";
    print $ofh "library(gplots)\n";
    print $ofh "library(Biobase)\n";
    
    
    print $ofh "load(\"$R_data_file\")\n";
    
    print $ofh "outdir = \"" . basename($R_data_file) . ".clusters_fixed_K_" . $K . "\"\n";
    print $ofh "dir.create(outdir)\n";
    
    print $ofh "gene_partition_assignments <- cutree(as.hclust(hc_genes), k=$K)\n";
    
    # make another heatmap:
    print $ofh "gene_colors = partition_colors[gene_partition_assignments]\n";
    print $ofh "postscript(file=\"clusters_fixed_K_${K}.heatmap.eps\", horizontal=FALSE, width=8, height=18, paper=\"special\");\n";
    print $ofh "heatmap.2(centered_data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors, scale=\"none\", density.info=\"none\", trace=\"none\", key=TRUE, keysize=1.2, cexCol=2.5, margins=c(15,15), lhei=c(0.3,2), lwid=c(2.5,4))\n";
    print $ofh "dev.off()\n";
    
    print $ofh "gene_names = rownames(centered_data)\n";
    print $ofh "num_cols = length(centered_data[1,])\n";
    
    
    print $ofh "for (i in 1:$K) {\n";
    print $ofh "    partition_i = (gene_partition_assignments == i)\n";
    
    print $ofh "    partition_centered_data = centered_data[partition_i,]\n";
    
    print $ofh "    # if the partition involves only one row, then it returns a vector instead of a table\n";
        ;
    print $ofh "    if (sum(partition_i) == 1) {\n";
    print $ofh "          dim(partition_centered_data) = c(1,num_cols)\n";
    print $ofh "          colnames(partition_centered_data) = colnames(centered_data)\n";
    print $ofh "          rownames(partition_centered_data) = gene_names[partition_i]\n";
    print $ofh "    }\n";
    
    
    print $ofh "    outfile = paste(outdir, \"/subcluster_\", i, sep='')\n";
    print $ofh "    write.table(partition_centered_data, file=outfile, quote=F, sep=\"\\t\")\n";
    print $ofh "}\n";
    

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
        die "Error, cmd $cmd died with ret $ret";
    }

    return;
}
