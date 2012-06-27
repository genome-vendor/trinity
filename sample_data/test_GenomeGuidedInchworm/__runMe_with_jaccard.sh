#!/bin/bash -ve

if [ -e bowtie_out.coordSorted.sam.gz ] && [ ! -e bowtie_out.coordSorted.sam ]; then
    gunzip -c bowtie_out.coordSorted.sam.gz > bowtie_out.coordSorted.sam
fi

## partition reads into non-overlapping coverage bins

../../util/prep_rnaseq_alignments_for_genome_assisted_assembly.pl --coord_sorted_SAM bowtie_out.coordSorted.sam -J 1 -I 500 --SS_lib_type RF --jaccard


## get list of filenames containing the partitioned reads

find Dir* -name "*reads" -type f > input.Read_files

## Create inchworm assembly commands based on the reads

../../util/write_inchworm_cmds.pl input.Read_files > iworm.cmds

## Run inchworm assemblies in parallel

../../Inchworm/bin/ParaFly -c iworm.cmds -CPU 4 -failed_cmds iworm.cmds.failed -v

## Combine all inchworm assemblies into a single output file, and assign unique accessions to each assembled contig.

find Dir* -type f -name "*iworm" -exec cat {} + | ../../util/inchworm_accession_incrementer.pl > iworm.GenomeGuided.fasta


## Done.

