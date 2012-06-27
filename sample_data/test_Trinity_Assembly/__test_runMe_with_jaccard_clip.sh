#!/bin/bash -ve

if [ -e reads.left.fq.gz ] && ! [ -e reads.left.fq ]
then
    gunzip -c reads.left.fq.gz > reads.left.fq
fi

if [ -e reads.right.fq.gz ] && ! [ -e reads.right.fq ]
then
    gunzip -c reads.right.fq.gz > reads.right.fq
fi


../../Trinity.pl --seqType fq --left reads.left.fq --right reads.right.fq --SS_lib_type RF  --CPU 4 --bfly_opts "-V 10 --stderr" --jaccard_clip --kmer_method meryl




