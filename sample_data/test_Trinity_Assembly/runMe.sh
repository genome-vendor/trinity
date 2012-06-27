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
../../Trinity.pl --seqType fq --kmer_method jellyfish --max_memory 2G --left reads.left.fq.gz --right reads.right.fq.gz --SS_lib_type RF --group_pairs_distance 501 --path_reinforcement_distance 76 --min_contig_length 201 --CPU 4 --bfly_opts "-V 10 --stderr" 

##### Done Running Trinity #####

if [ ! $* ]; then
    exit 0
fi


sleep 2

######################################################
## align reads back to the transcripts using Bowtie ##
######################################################

sleep 2

../../util/alignReads.pl --left reads.left.fq --right reads.right.fq --target trinity_out_dir/Trinity.fasta --aligner bowtie --seqType fq --SS_lib_type RF

##### Done aligning reads #######

sleep 2

###########################################
# use RSEM to estimate read abundance  ####
###########################################

sleep 2

../../util/RSEM_util/run_RSEM.pl --transcripts trinity_out_dir/Trinity.fasta --name_sorted_bam bowtie_out/bowtie_out.nameSorted.sam.+.sam.PropMapPairsForRSEM.bam --paired 

###### Done running RSEM ########

sleep 2


############################
# compute FPKM values   ####
############################

sleep 2

../../util/RSEM_util/summarize_RSEM_fpkm.pl --transcripts trinity_out_dir/Trinity.fasta --RSEM RSEM.isoforms.results --fragment_length 300 --group_by_component | tee Trinity.RSEM.fpkm


#############################
####   Done.  ###############
#############################




