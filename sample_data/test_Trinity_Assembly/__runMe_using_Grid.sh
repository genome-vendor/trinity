#!/bin/bash -ve

if [ -e reads.left.fq.gz ] && ! [ -e reads.left.fq ]
then
    gunzip -c reads.left.fq.gz > reads.left.fq
fi

if [ -e reads.right.fq.gz ] && ! [ -e reads.right.fq ]
then
    gunzip -c reads.right.fq.gz > reads.right.fq
fi


#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################



## use jellyfish
../../Trinity.pl --seqType fq --kmer_method jellyfish --max_memory 2G --left reads.left.fq.gz --right reads.right.fq.gz --SS_lib_type RF --paired_fragment_length 300  --min_contig_length 200 --CPU 4 --bfly_opts "-V 10 --stderr" --grid_computing_module BroadInstGridRunner



