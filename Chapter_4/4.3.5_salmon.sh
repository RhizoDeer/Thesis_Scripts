#!/bin/bash

#SBATCH -o slurm-%x_%A_%a.out 
#SBATCH --cpus-per-task=8
#SBATCH --mem=140G
#SBATCH --partition=medium
#SBATCH --array=1-24


##########################################################################
#@author Micha Bayer, James Hutton Institute
#Modified Katie Arnton, University of Dundee, Mar 2026
# Configure and start
##########################################################################

#absolute path to the folder containing the index folder for Salmon
salmonIndexPath=/mnt/shared/projects/jhi/barley/202603_RNASeqJH203/BarkeRTD
#the name of the folder within the above that contains the actual index files -- we will need to copy this to the local scratch disk
salmonIndexFolder=BarkeRTD_salmon_index
#a list of paths to the R1 files of all samples
sampleList=allR1Files.txt

echo "starting run"
date

#list the host and tmp dir
echo "host = $HOSTNAME"
echo -e "\nTMPDIR = $TMPDIR"

##########################################################################
# Pick a sample from the list
##########################################################################

#parse an input file with the sample paths
#split the input by line
IFS=$'\n'
#read the whole file 
lines=( `cat "$sampleList" `)
#based on the SLURM_ARRAY_TASK_ID we pick a single line from the input file, and this will be the line that gets used
#note the SLURM_ARRAY_TASK_ID starts at 1, but the array is zero-based, hence need to subtract 1 
R1FilePath=${lines[$SLURM_ARRAY_TASK_ID-1]}
#derive the name of the R2 file
R2FilePath=`echo $R1FilePath | sed 's/R1_001.fastq.gz/R2_001.fastq.gz/g'`

#work out the sample name from the path
sampleName=`basename $R1FilePath _R1_001.fastq.gz`

echo -e "\n\n============================"
echo -e "processing sample $sampleName"
echo "R1FilePath = $R1FilePath"
echo "R2FilePath = $R2FilePath"

#make a dedicated dir for this
mkdir $sampleName
cd $sampleName


###########################################################################
# Quality/adapter trim the reads
###########################################################################
echo -e "\nrun fastp"
date

#fastp has its own conda environment
source activate fastp

fastp \
-i $R1FilePath \
-I $R2FilePath \
-o $sampleName.out.R1.fq \
-O $sampleName.out.R2.fq \
--html $sampleName.html \
--json $sampleName.json \
-q 20 \
--cut_front \
--cut_tail \
--trim_poly_g \
-l 30 \
--detect_adapter_for_pe \
--thread $SLURM_CPUS_PER_TASK \
>> $sampleName.fastp.log 2>&1

conda deactivate

###########################################################################
# Quantify the reads with Salmon
###########################################################################
source activate OldSalmon

echo -e "\nrun Salmon"
date
salmon quant \
-i $salmonIndexPath/$salmonIndexFolder \
-l A \
-1 $sampleName.out.R1.fq \
-2 $sampleName.out.R2.fq \
-p $SLURM_CPUS_PER_TASK \
--seqBias \
--gcBias \
--posBias \
--validateMappings \
-o . \
>> $sampleName.salmon.log 2>&1


echo "workflow complete"
date






