#!/bin/bash -ve

if [ -e reads.left.fq.gz ] && ! [ -e reads.left.fq ]
then
    gunzip -c reads.left.fq.gz > reads.left.fq
fi


../../Trinity.pl --seqType fq --single reads.left.fq  --min_contig_length 305 --kmer_method meryl  --CPU 4 --bfly_opts "-V 10 --stderr"  --output trinity_single_outdir




