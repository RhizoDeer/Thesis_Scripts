#!/bin/bash

#SBATCH -o slurm-%x_%A.out
#SBATCH --mem=150G
#SBATCH --cpus-per-task=12

############################################################################################
#Script for creating decoy-aware Salmon index files
#We will be following the approach that uses entire genomes as decoys as this is more accurate
#Follows instructions at https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

#@author Micha Bayer, James Hutton Institute, Sep 2021
#Modified Katie Arnton, University of Dundee, Mar 2026
############################################################################################


#=============================
#VARIABLES
#=============================

#a FASTA file containing the reference transcriptome
transcriptomeFASTA=/mnt/shared/projects/jhi/barley/202603_RNASeqJH203/BaRT2v18.fa
#a FASTA file containg the corresponding reference genome
genomeFASTA=/mnt/shared/projects/jhi/barley/202603_RNASeqJH203/180903_Barke_Unfiltered_chloro_clean_pseudomolecules_v1.fasta
#a prefix for labeling our output
prefix=BarkeRTD
#the number of threads we want to use
numThreads=$SLURM_CPUS_PER_TASK

#do everything in a dedicated directory
mkdir $prefix
cd $prefix

#extract the names of the reference sequences in the genome file
echo "make $prefix.decoys.txt file"
grep "^>" $genomeFASTA \
| cut -d " " -f 1 \
> $prefix.decoys.txt

#get rid of the ">" characters
sed -i.bak -e 's/>//g' $prefix.decoys.txt

#concatenate the transcriptome and genome FASTA files, in that order (important)
echo "make combined FASTA file"
cat \
$transcriptomeFASTA \
$genomeFASTA \
> $prefix.withDecoys.fasta

#Activate conda environment w/ Salmon version 1.10.2
source activate OldSalmon

#now we are ready to build the Salmon index
echo "make Salmon index"
salmon index \
-t $prefix.withDecoys.fasta \
-d $prefix.decoys.txt \
-p $numThreads \
-i $prefix"_salmon_index"

echo "workflow complete"
