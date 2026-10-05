#####################################################################################
#HYPOTHESIS TESTING part II
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
library("UpSetR")

#############################################################
#set working directory-Davide CPU
setwd("/cluster/db/R_shared/JH201/")
#set working directory-Katie CPU
setwd("/cluster/db/carnton/JH201")
getwd()
#############################################################

##################################################################
#Import the pre-processed DESeq2 objects: RKMP
#################################################################

#import the datasets
#JH201_RKMP_Tube_DESeq <-readRDS("JH201_RKMP_tube.rds")
JH201_RKMP_Tube_DESeq <-readRDS("Original_data/JH201_RKMP_tube.rds")


#inspect the data
JH201_RKMP_Tube_DESeq

#execute the differential count analysis with the function DESeq 
JH201_RKMP_Tube_test <- DESeq(JH201_RKMP_Tube_DESeq, fitType="local", betaPrior = FALSE) 

#################################
#Plant effect: genotypes vs. Bulk
#################################

#Barke
Barke_vs_unplanted <- results(JH201_RKMP_Tube_test, contrast = c("Genotype", "Barke", "Bulk")) 

#Inspect the files
Barke_vs_unplanted[1:6, ]
#these files contains the output of the analysis (indeed the results), we need to filter for the one below the significant threshold imposed and enriched in the plant microhabitats
#operationally this is a subsetting of the original results files

#extract  ASVs whose adjusted p.value in a given comparison is below 0.01
#what is an 'adjusted p value' and what kind of adjustment we have to make? https://en.wikipedia.org/wiki/Multiple_comparisons_problem 
Barke_vs_unplanted_FDR001 <- Barke_vs_unplanted[(rownames(Barke_vs_unplanted)[which(Barke_vs_unplanted$padj <0.01)]), ]
Barke_vs_unplanted_FDR001

#enriched in Barke: Barke is the first term of comparison, positive fold changes
Barke_enriched <-  Barke_vs_unplanted_FDR001[(rownames(Barke_vs_unplanted_FDR001)[which(Barke_vs_unplanted_FDR001$log2FoldChange > 0)]), ]
Barke_enriched

#enriched in Bulk: Barke is the first term of comparison, negative fold changes
Bulk_enriched <-  Barke_vs_unplanted_FDR001[(rownames(Barke_vs_unplanted_FDR001)[which(Barke_vs_unplanted_FDR001$log2FoldChange < 0)]), ]
Bulk_enriched

#Morex
Morex_vs_unplanted <- results(JH201_RKMP_Tube_test, contrast = c("Genotype", "Morex", "Bulk")) 

#Inspect the files
Morex_vs_unplanted[1:6, ]

#extract  ASVs whose adjusted p.value in a given comparison is below 0.01
Morex_vs_unplanted_FDR001 <- Morex_vs_unplanted[(rownames(Morex_vs_unplanted)[which(Morex_vs_unplanted$padj <0.01)]), ]
Morex_vs_unplanted_FDR001

#enriched in Morex: Morex is the first term of comparison, positive fold changes
Morex_enriched <-  Morex_vs_unplanted_FDR001[(rownames(Morex_vs_unplanted_FDR001)[which(Morex_vs_unplanted_FDR001$log2FoldChange > 0)]), ]
Morex_enriched

#enriched in Bulk: Morex is the first term of comparison, negative fold changes
Bulk_enriched <-  Morex_vs_unplanted_FDR001[(rownames(Morex_vs_unplanted_FDR001)[which(Morex_vs_unplanted_FDR001$log2FoldChange < 0)]), ]
Bulk_enriched

#SL17
SL17_vs_unplanted <- results(JH201_RKMP_Tube_test, contrast = c("Genotype", "SL17", "Bulk")) 

#Inspect the files
SL17_vs_unplanted[1:6, ]

#extract  ASVs whose adjusted p.value in a given comparison is below 0.01
SL17_vs_unplanted_FDR001 <- SL17_vs_unplanted[(rownames(SL17_vs_unplanted)[which(SL17_vs_unplanted$padj <0.01)]), ]
SL17_vs_unplanted_FDR001

#enriched in SL17: SL17 is the first term of comparison, positive fold changes
SL17_enriched <-  SL17_vs_unplanted_FDR001[(rownames(SL17_vs_unplanted_FDR001)[which(SL17_vs_unplanted_FDR001$log2FoldChange > 0)]), ]
SL17_enriched

#enriched in Bulk: SL17 is the first term of comparison, negative fold changes
Bulk_enriched <-  SL17_vs_unplanted_FDR001[(rownames(SL17_vs_unplanted_FDR001)[which(SL17_vs_unplanted_FDR001$log2FoldChange < 0)]), ]
Bulk_enriched

#SL52
SL52_vs_unplanted <- results(JH201_RKMP_Tube_test, contrast = c("Genotype", "SL52", "Bulk")) 

#Inspect the files
SL52_vs_unplanted[1:6, ]

#extract  ASVs whose adjusted p.value in a given comparison is below 0.01
SL52_vs_unplanted_FDR001 <- SL52_vs_unplanted[(rownames(SL52_vs_unplanted)[which(SL52_vs_unplanted$padj <0.01)]), ]
SL52_vs_unplanted_FDR001

#enriched in SL52: SL52 is the first term of comparison, positive fold changes
SL52_enriched <-  SL52_vs_unplanted_FDR001[(rownames(SL52_vs_unplanted_FDR001)[which(SL52_vs_unplanted_FDR001$log2FoldChange > 0)]), ]
SL52_enriched

#enriched in Bulk: SL52 is the first term of comparison, negative fold changes
Bulk_enriched <-  SL52_vs_unplanted_FDR001[(rownames(SL52_vs_unplanted_FDR001)[which(SL52_vs_unplanted_FDR001$log2FoldChange < 0)]), ]
Bulk_enriched

##########################################
#Genotype effect: comparisons against Barke
##########################################

#Barke vs SL52
Barke_vs_SL52 <- results(JH201_RKMP_Tube_test, contrast = c("Genotype", "Barke", "SL52")) 

#Inspect the files
Barke_vs_SL52[1:6, ]

#extract  ASVs whose adjusted p.value in a given comparison is below 0.01
Barke_vs_SL52_FDR001 <- Barke_vs_SL52[(rownames(Barke_vs_SL52)[which(Barke_vs_SL52$padj <0.01)]), ]
Barke_vs_SL52_FDR001

#Barke vs SL17
Barke_vs_SL17 <- results(JH201_RKMP_Tube_test, contrast = c("Genotype", "Barke", "SL17")) 

#Inspect the files
Barke_vs_SL17[1:6, ]

#extract  ASVs whose adjusted p.value in a given comparison is below 0.01
Barke_vs_SL17_FDR001 <- Barke_vs_SL17[(rownames(Barke_vs_SL17)[which(Barke_vs_SL17$padj <0.01)]), ]
Barke_vs_SL17_FDR001

#Barke vs Morex
Barke_vs_Morex <- results(JH201_RKMP_Tube_test, contrast = c("Genotype", "Barke", "Morex")) 

#Inspect the files
Barke_vs_Morex[1:6, ]

#extract  ASVs whose adjusted p.value in a given comparison is below 0.01
Barke_vs_Morex_FDR001 <- Barke_vs_Morex[(rownames(Barke_vs_Morex)[which(Barke_vs_Morex$padj <0.01)]), ]
Barke_vs_Morex_FDR001

#####################################################################
#Conclusions
# 1) the Plant effect observed in ordination methods appears sustained by the enrichment of Pedobacter
# 2) Morex emerged with a distinct profile as it enriches another taxa, haven't checked but it is likely to be Pseudomonas
# 3) No individual enrichment differentiate between genotypes, so what we have seen in PERMANOVA is a combination of bacteria, not just one
#####################################################################

##################################################################
#Import the pre-processed DESeq2 objects: no27
#################################################################

#import the datasets
JH201_no27_Tube_DESeq <-readRDS("Original_data/JH201_no27_tube.rds")

#inspect the data
JH201_no27_Tube_DESeq

#execute the differential count analysis with the function DESeq 
JH201_no27_Tube_test <- DESeq(JH201_no27_Tube_DESeq, fitType="local", betaPrior = FALSE) 

####################
#Plant effect: genotypes vs. Bulk
####################

#Barke
Barke_vs_unplanted <- results(JH201_no27_Tube_test, contrast = c("Genotype", "Barke", "Bulk")) 

#Inspect the files
Barke_vs_unplanted[1:6, ]

#extract  ASVs whose adjusted p.value in a given comparison is below 0.01
Barke_vs_unplanted_FDR001 <- Barke_vs_unplanted[(rownames(Barke_vs_unplanted)[which(Barke_vs_unplanted$padj <0.01)]), ]
Barke_vs_unplanted_FDR001

#enriched in Barke: Barke is the first term of comparison, positive fold changes
Barke_enriched <-  Barke_vs_unplanted_FDR001[(rownames(Barke_vs_unplanted_FDR001)[which(Barke_vs_unplanted_FDR001$log2FoldChange > 0)]), ]
Barke_enriched

#enriched in Bulk: Barke is the first term of comparison, negative fold changes
Bulk_enriched <-  Barke_vs_unplanted_FDR001[(rownames(Barke_vs_unplanted_FDR001)[which(Barke_vs_unplanted_FDR001$log2FoldChange < 0)]), ]
Bulk_enriched

#Morex
Morex_vs_unplanted <- results(JH201_no27_Tube_test, contrast = c("Genotype", "Morex", "Bulk")) 

#Inspect the files
Morex_vs_unplanted[1:6, ]

#extract  ASVs whose adjusted p.value in a given comparison is below 0.01
Morex_vs_unplanted_FDR001 <- Morex_vs_unplanted[(rownames(Morex_vs_unplanted)[which(Morex_vs_unplanted$padj <0.01)]), ]
Morex_vs_unplanted_FDR001

#enriched in Morex: Morex is the first term of comparison, positive fold changes
Morex_enriched <-  Morex_vs_unplanted_FDR001[(rownames(Morex_vs_unplanted_FDR001)[which(Morex_vs_unplanted_FDR001$log2FoldChange > 0)]), ]
Morex_enriched

#enriched in Bulk: Morex is the first term of comparison, negative fold changes
Bulk_enriched <-  Morex_vs_unplanted_FDR001[(rownames(Morex_vs_unplanted_FDR001)[which(Morex_vs_unplanted_FDR001$log2FoldChange < 0)]), ]
Bulk_enriched

#SL17
SL17_vs_unplanted <- results(JH201_no27_Tube_test, contrast = c("Genotype", "SL17", "Bulk")) 

#Inspect the files
SL17_vs_unplanted[1:6, ]

#extract  ASVs whose adjusted p.value in a given comparison is below 0.01
SL17_vs_unplanted_FDR001 <- SL17_vs_unplanted[(rownames(SL17_vs_unplanted)[which(SL17_vs_unplanted$padj <0.01)]), ]
SL17_vs_unplanted_FDR001

#enriched in SL17: SL17 is the first term of comparison, positive fold changes
SL17_enriched <-  SL17_vs_unplanted_FDR001[(rownames(SL17_vs_unplanted_FDR001)[which(SL17_vs_unplanted_FDR001$log2FoldChange > 0)]), ]
SL17_enriched

#enriched in Bulk: SL17 is the first term of comparison, negative fold changes
Bulk_enriched <-  SL17_vs_unplanted_FDR001[(rownames(SL17_vs_unplanted_FDR001)[which(SL17_vs_unplanted_FDR001$log2FoldChange < 0)]), ]
Bulk_enriched

#SL52
SL52_vs_unplanted <- results(JH201_no27_Tube_test, contrast = c("Genotype", "SL52", "Bulk")) 

#Inspect the files
SL52_vs_unplanted[1:6, ]

#extract  ASVs whose adjusted p.value in a given comparison is below 0.01
SL52_vs_unplanted_FDR001 <- SL52_vs_unplanted[(rownames(SL52_vs_unplanted)[which(SL52_vs_unplanted$padj <0.01)]), ]
SL52_vs_unplanted_FDR001

#enriched in SL52: SL52 is the first term of comparison, positive fold changes
SL52_enriched <-  SL52_vs_unplanted_FDR001[(rownames(SL52_vs_unplanted_FDR001)[which(SL52_vs_unplanted_FDR001$log2FoldChange > 0)]), ]
SL52_enriched

#enriched in Bulk: SL52 is the first term of comparison, negative fold changes
Bulk_enriched <-  SL52_vs_unplanted_FDR001[(rownames(SL52_vs_unplanted_FDR001)[which(SL52_vs_unplanted_FDR001$log2FoldChange < 0)]), ]
Bulk_enriched

##########################################
#Genotype effect: comparisons against Barke
##########################################

#Barke vs SL52
Barke_vs_SL52 <- results(JH201_no27_Tube_test, contrast = c("Genotype", "Barke", "SL52")) 

#Inspect the files
Barke_vs_SL52[1:6, ]

#extract  ASVs whose adjusted p.value in a given comparison is below 0.01
Barke_vs_SL52_FDR001 <- Barke_vs_SL52[(rownames(Barke_vs_SL52)[which(Barke_vs_SL52$padj <0.01)]), ]
Barke_vs_SL52_FDR001

#Barke vs SL17
Barke_vs_SL17 <- results(JH201_no27_Tube_test, contrast = c("Genotype", "Barke", "SL17")) 

#Inspect the files
Barke_vs_SL17[1:6, ]

#extract  ASVs whose adjusted p.value in a given comparison is below 0.01
Barke_vs_SL17_FDR001 <- Barke_vs_SL17[(rownames(Barke_vs_SL17)[which(Barke_vs_SL17$padj <0.01)]), ]
Barke_vs_SL17_FDR001

#Barke vs Morex
Barke_vs_Morex <- results(JH201_no27_Tube_test, contrast = c("Genotype", "Barke", "Morex")) 

#Inspect the files
Barke_vs_Morex[1:6, ]

#extract  ASVs whose adjusted p.value in a given comparison is below 0.01
Barke_vs_Morex_FDR001 <- Barke_vs_Morex[(rownames(Barke_vs_Morex)[which(Barke_vs_Morex$padj <0.01)]), ]
Barke_vs_Morex_FDR001

##################################################################
#Import the pre-processed DESeq2 objects: Barke & RKMP
#################################################################

#import the datasets
JH201_RKMP_Barke_Tube_DESeq <-readRDS("Original_data/JH201_RKMP_Barke_tube.rds")

#inspect the data
JH201_RKMP_Barke_Tube_DESeq

#execute the differential count analysis with the function DESeq 
JH201_RKMP_Barke_Tube_test <- DESeq(JH201_RKMP_Barke_Tube_DESeq, fitType="local", betaPrior = FALSE)

####################
#Run effect: 230406 vs. 230525 and 230831
####################

#230525
r230406_vs_r230525 <- results(JH201_RKMP_Barke_Tube_test, contrast = c("Date.Inoc", "230406", "230525")) 

#Inspect the files
r230406_vs_r230525[1:6, ]

#extract  ASVs whose adjusted p.value in a given comparison is below 0.01
r230406_vs_r230525_FDR001 <- r230406_vs_r230525[(rownames(r230406_vs_r230525)[which(r230406_vs_r230525$padj <0.01)]), ]
r230406_vs_r230525_FDR001

#enriched in 230406: 230406 is the first term of comparison, positive fold changes
r230406_enriched <-  r230406_vs_r230525_FDR001[(rownames(r230406_vs_r230525_FDR001)[which(r230406_vs_r230525_FDR001$log2FoldChange > 0)]), ]
r230406_enriched

#enriched in 230525: 230406 is the first term of comparison, negative fold changes
r230525_enriched <-  r230406_vs_r230525_FDR001[(rownames(r230406_vs_r230525_FDR001)[which(r230406_vs_r230525_FDR001$log2FoldChange < 0)]), ]
r230525_enriched


#230831
r230406_vs_r230831 <- results(JH201_RKMP_Barke_Tube_test, contrast = c("Date.Inoc", "230406", "230831")) 

#Inspect the files
r230406_vs_r230831[1:6, ]

#extract  ASVs whose adjusted p.value in a given comparison is below 0.01
r230406_vs_r230831_FDR001 <- r230406_vs_r230831[(rownames(r230406_vs_r230831)[which(r230406_vs_r230831$padj <0.01)]), ]
r230406_vs_r230831_FDR001

#####################################################################
#double check sample distribution and counts
#####################################################################

