#####################################################################################
#PREPROCESING
####################################################################################
#############################################################
#
# Ref to the ARTICLE 
# 
#  Code to compute calculations presented in Katie's work
#  Revision 02/24
#  carnton@dundee.ac.uk
#  d.bulgarelli@dundee.ac.uk 

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################
#required packages 
library("phyloseq")

#############################################################
#set working directory: DB cpu
#setwd("/cluster/db/R_shared/JH201/")
#set working directory: Katie
setwd("/cluster/db/carnton/JH201")
#check working directory
getwd()

#import the dataset
#JH201 <-readRDS("JH201_dada2_silva_138.1.rds")
JH201 <-readRDS("Original_data/JH201_dada2_silva_138.1.rds")

#inspect the files
JH201

#and its three "constituents"

#ASV counts
otu_table(JH201)

#Taxonomy information
tax_table(JH201)

#mapping files
sample_data(JH201)

#check the sample names
sample_names(JH201)

#export the files to simplify the mapping file
#write.table(as.data.frame(sample_data(JH201)), file = "JH201_map1.txt")

#import the simplified mapping files; note the use of the "corrected" version
#JH201_new_map <- read.delim("JH201_map2.txt", row.names = 1)
JH201_new_map <- read.delim("Original_data/JH201_map2.txt", row.names = 1)
JH201_new_map

#and replace the old mapping file
sample_data(JH201) <- JH201_new_map
sample_data(JH201)

#Example of subsetting: positive selection
JH201_inoculant <- subset_samples(JH201,  Sample.Type == "Inoculant")
JH201_inoculant
#Example of subsetting: negative election (when more than two levels for a factor)
JH201_noPBS <- subset_samples(JH201,  SynCom.ID != "PBS")
JH201_noPBS
sample_data(JH201_noPBS)

#stats overall library
sum(sample_sums(JH201))

#proceed with independent filtering in both objects prior merging
##################################################################
#Pre-processing: remove Chloroplast and Mitochondria but retain NA
#rationale: 16S rRNA primers may amplify plant-derived sequences 
#those may interfere with data analysis (we are after bacteria, not plants)
#################################################################

#JH201
JH201_no_chlor <-subset_taxa(JH201, (Order!="Chloroplast") | is.na(Order))
JH201_no_chlor

JH201_no_plants <-subset_taxa(JH201_no_chlor, (Family!="Mitochondria") | is.na(Family))
JH201_no_plants

#########################################################################
# Remove ASVs assigned to NA at phylum level
# Rationale: if a sequence cannot be assigned at a very high level such as Phylum is not of use us
#########################################################################

JH201_no_plants_1 <- subset_taxa(JH201_no_plants, Phylum!= "NA")
JH201_no_plants_1

#########################################################################
# inspect reads distribution across objects
# Rationale: if we want to compare like-with-like we should know sample distribution
#########################################################################

#summary
sum (sample_sums(JH201_no_plants_1))
min(sample_sums(JH201_no_plants_1))
max(sample_sums(JH201_no_plants_1))
mean(sample_sums(JH201_no_plants_1))

#graphical outputs
sort(sample_sums(JH201_no_plants_1))
hist(sample_sums(JH201_no_plants_1), main = paste("JH201 reads distribution"), xlab = paste("reads"), ylab = paste ("number of samples"))

##############################################################################################################################################
# Filter low abundance ASVs *TO BE CHECKED BY DB*
# ASVs represented by few reads in one or few samples may represents either the so-called rare biosphere or artefacts of the sequencing approach
#Either way they are poorly reproducible and confound the statistical analyses
#We operate a so called 'secondary filtering approach' whose rationale is described here https://www.nature.com/articles/nmeth.2276 
#############################################################################################################################################

#Set an arbitrary threshold of 10 reads in 20% of samples
JH201_threshold <- filter_taxa(JH201_no_plants_1, function(x) sum(x > 10) > (0.2*length(x)), TRUE)
#example to tweak the parameters, e.g., 5 reads in 10% of samples
#JH201_threshold <- filter_taxa(JH201_no_plants_1, function(x) sum(x > 5) > (0.1*length(x)), TRUE)

#The impact of secondary filtering
JH201_no_plants_1
JH201_threshold

#wow, that's a massive filtering: we retained little more than 6% of the initial ASVs
#measure the impact of thresholding
sort(sample_sums(JH201_threshold))
hist(sample_sums(JH201_threshold))

##ratio filtered reads/total reads
ratio <- sum(sample_sums(JH201_threshold))/sum(sample_sums(JH201_no_plants_1))*100
ratio

#summary
min(sample_sums(JH201_threshold))
max(sample_sums(JH201_threshold))
mean(sample_sums(JH201_threshold))

###################################################################################################
#agglomerate samples at genus level
###################################################################################################

JH201_genus_10_reads <- tax_glom(JH201_threshold, taxrank='Genus')
JH201_genus_10_reads
#Inspect the genera
tax_table(JH201_genus_10_reads)

###################################################################################################
#Generate .rds objects for downstream analyses
#more info about .rds format https://riptutorial.com/r/example/3650/rds-and-rdata--rda--files
###################################################################################################

#saveRDS(JH201_genus_10_reads, file = "JH201_genus_10_reads.rds")
###################################################################################################
#Key points
###################################################################################################

#1) The three constituent files of a Phyloseq object
#2) The rationale of pre-processing: a sequencing data "as is" may not have biological value
#3) the .rds format(though this concept will become lot clearer in the next session)

###################################################################################################

#Adjust to lower threshold

#50 reads in 5% of the samples
JH201_threshold_50r5p <- filter_taxa(JH201_no_plants_1, function(x) sum(x > 50) > (0.05*length(x)), TRUE)
JH201_threshold_50r5p
ratio <- sum(sample_sums(JH201_threshold_50r5p))/sum(sample_sums(JH201_no_plants_1))*100
ratio
JH201_genus_50r5p <- tax_glom(JH201_threshold_50r5p, taxrank='Genus')
JH201_genus_50r5p
tax_table(JH201_genus_50r5p)

#saveRDS(JH201_genus_50r5p, file = "Original_data/JH201_genus_50_reads.rds")

####Inoculum alone####

#Subset prior to thresholds for the low abundance ASVs
JH201_inoculant <- subset_samples(JH201_no_plants_1,  Sample.Type == "Inoculant")
JH201_inoculant

#view sample variables
sample_data(JH201_inoculant)

#Set threshold for ASV filtering at 20% and 10 reads
JH201_inoculant_threshold_10r20p <- filter_taxa(JH201_inoculant, function(x) sum(x > 10) > (0.2*length(x)), TRUE)

#Compare with original
JH201_inoculant
JH201_inoculant_threshold_10r20p

#measure the impact of thresholding
sort(sample_sums(JH201_inoculant_threshold_10r20p))
hist(sample_sums(JH201_inoculant_threshold_10r20p))

##ratio filtered reads/total reads
ratio <- sum(sample_sums(JH201_inoculant_threshold_10r20p))/sum(sample_sums(JH201_inoculant))*100
ratio

#summary
min(sample_sums(JH201_inoculant_threshold_10r20p))
max(sample_sums(JH201_inoculant_threshold_10r20p))
mean(sample_sums(JH201_inoculant_threshold_10r20p))

#agglomerate
JH201_inoculant_genus_10r20p <- tax_glom(JH201_inoculant_threshold_10r20p, taxrank='Genus')
JH201_inoculant_genus_10r20p

#Inspect the genera
tax_table(JH201_inoculant_genus_10r20p)

#Save processed reads for inoculum
#saveRDS(JH201_inoculant_genus_10r20p, file = "Original_data/JH201_inoculant_genus_10r20p.rds")
