#!/bin/sh
#PBS -v PATH
#$ -v PATH


para=$1
cd /home/aoli/Documents/RBP_surface_make_database/data/cd-hit_0906/03
./RBP04.fasta.50319-bl.pl 0 $para &
wait

