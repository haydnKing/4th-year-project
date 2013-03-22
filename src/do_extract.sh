#!/bin/bash
for file in Genomes/*.fasta
do
    python extract.py $file
done
