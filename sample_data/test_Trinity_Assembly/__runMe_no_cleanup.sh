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
../../Trinity.pl --seqType fq --kmer_method jellyfish --max_memory 2G --left reads.left.fq.gz --right reads.right.fq.gz --SS_lib_type RF --group_pairs_distance 501 --path_reinforcement_distance 76 --min_contig_length 201 --CPU 4 --no_cleanup --bfly_opts "-V 10 --stderr" 

##### Done Running Trinity #####


exit 0
