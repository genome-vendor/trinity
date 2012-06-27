#!/bin/bash -v

example_100trans.coord_sorted.sam.gz  example_100trans.fa.gz  example_100trans.name_sorted.sam

if [ -e example_100trans.fa.gz ] && ! [ -e example_100trans.fa ]
then
    gunzip -c example_100trans.fa.gz > example_100trans.fa
fi

if [ -e example_100trans.coord_sorted.sam.gz ] && ! [ -e example_100trans.coord_sorted.sam ]
then
    gunzip -c example_100trans.coord_sorted.sam.gz > example_100trans.coord_sorted.sam
fi

if ! [ -e example_100trans.name_sorted.sam ]
then
    sort -T . -S 2G -k 1,1 -k 3,3 example_100trans.coord_sorted.sam > example_100trans.name_sorted.sam
fi


# run Trinity_EM

../../Inchworm/bin/Trinity_EM --fasta example_100trans.fa --name_sorted_sam example_100trans.name_sorted.sam --coord_sorted_sam example_100trans.coord_sorted.sam --min_delta_ML 0.0001 --fragment_length 1


