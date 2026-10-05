#####################################################################################
#HYPOTHESIS TESTING
####################################################################################
#############################################################
#
# Ref to the ARTICLE 
# 
#  Code to compute calculations presented in Katie's work
#  Revision 01/24
#  c.arnton@dundee.ac.uk
#  d.bulgarelli@dundee.ac.uk 

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################
library("phyloseq")
library("DESeq2")

#############################################################
#set working directory-Davide CPU
#setwd("/cluster/db/R_shared/JH201/")
#set working directory-Katie CPU
setwd("/cluster/db/carnton/JH201")
getwd()
#############################################################

##################################################################
#Import the pre-processed Phyloseq object
#################################################################

#import the dataset: integer counts
#JH201 <-readRDS("JH201_genus_50_reads.rds")
JH201 <-readRDS("Original_data/JH201_genus_50_reads.rds")
#more info about .rds format https://riptutorial.com/r/example/3650/rds-and-rdata--rda--files

#inspect the files
JH201
sample_data(JH201)

###################################################################################################
#Hypothesis testing
# Rationale: in the previous code we identified a significant effect of the two main independent variables we are looking at: SynCom and Genotype
# Here we ask the question how many ASVs (and what taxa) underpin this diversification 
# This will be critical to categorise differences between independent variables.
# The underlying assumption in our calculation is that, for being classified as a 'plant effect', we require for a given feature (i.e., a genus) to be signficantly enriched compared to Bulk
# Due to a) the distribution of the data and b) the experimental design we will adopt an approach called DESeq
# Here is additional  info: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y
# Here is the reference to the 'Sietske van Bentum's street' presented at Symposium https://www.nature.com/articles/s41467-022-28034-z
# Further information on method choice can be also found in the rebuttal letter of Carmen's paper
# This code is comprised by two parts: an object preparation and the actual calculation/data visualisation
##################################################################################################

#split the phyloseq object: starting to run with full set of genotypes
JH201_RKMP <- subset_samples(JH201, Date.Inoc  == "230406")
JH201_no27 <- subset_samples(JH201, Date.Inoc == "230223")
JH201_both <- subset_samples(JH201, Date.Inoc == "230525")

####################################
#RKMP
###################################

#prune non tube samples
JH201_RKMP_tube <- subset_samples(JH201_RKMP, Sample.Type  == "Tube")
sample_data(JH201_RKMP_tube)

#prepare the DESeq2 object: it requires two files called count and col data which will be extracted from the Phyloseq object

#extract count data 
JH201_RKMP_tube_counts <- otu_table(JH201_RKMP_tube)
countData = as.data.frame(JH201_RKMP_tube_counts)
colnames(JH201_RKMP_tube_counts)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(JH201_RKMP_tube)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
JH201_RKMP_tube_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)
#note on setting the formula
#saveRDS(JH201_RKMP_tube_cds, file = "Original_data/JH201_RKMP_tube.rds")

####################################
#no Pedobacter SynCom
###################################

#prune non tube samples
JH201_no27_tube <- subset_samples(JH201_no27, Sample.Type  == "Tube")
sample_data(JH201_no27_tube)

#prepare the DESeq2 object

#extract count data 
JH201_no27_tube_counts <- otu_table(JH201_no27_tube)
countData = as.data.frame(JH201_no27_tube_counts)
colnames(JH201_no27_tube_counts)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(JH201_no27_tube)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
JH201_no27_tube_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)
#note on setting the formula
#saveRDS(JH201_no27_tube_cds, file = "Original_data/JH201_no27_tube.rds")

####################################
#Both SynComs
###################################

#prune non tube samples and controls
JH201_both_tube <- subset_samples(JH201_both, Sample.Type  == "Tube")
JH201_both_tube <- subset_samples(JH201_both_tube, SynCom.ID  != "SDW")
sample_data(JH201_both_tube)

#prepare the DESeq2 object

#extract count data 
JH201_both_tube_counts <- otu_table(JH201_both_tube)
countData = as.data.frame(JH201_both_tube_counts)
colnames(JH201_both_tube_counts)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(JH201_both_tube)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
JH201_both_tube_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)
#note on setting the formula
#saveRDS(JH201_both_tube_cds, file = "Original_data/JH201_both_tube.rds")

####################################
#Barke and RKMP
###################################

#prune non tube samples and controls
JH201_RKMP_Barke_tube <- subset_samples(JH201, Sample.Type  == "Tube")
JH201_RKMP_Barke_tube <- subset_samples(JH201_RKMP_Barke_tube, SynCom.ID  == "RKMP")
JH201_RKMP_Barke_tube <- subset_samples(JH201_RKMP_Barke_tube, Genotype  == "Barke")
sample_data(JH201_RKMP_Barke_tube)

#prepare the DESeq2 object

#extract count data 
JH201_RKMP_Barke_tube_counts <- otu_table(JH201_RKMP_Barke_tube)
countData = as.data.frame(JH201_RKMP_Barke_tube_counts)
colnames(JH201_RKMP_Barke_tube_counts)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(JH201_RKMP_Barke_tube)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
JH201_RKMP_Barke_tube_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Date.Inoc)
#note on setting the formula
#saveRDS(JH201_RKMP_Barke_tube_cds, file = "Original_data/JH201_RKMP_Barke_tube.rds")
