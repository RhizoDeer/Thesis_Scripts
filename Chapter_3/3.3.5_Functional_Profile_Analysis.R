#############################################################
#############################################################
#
#  Analysis for chapter 3
#  Revision 08/26
#  c.arnton@dundee.ac.uk
#
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
library("vegan")
library ("ggplot2")
library ("DESeq2")
library("viridis")
library("gridExtra")
library("readr")
library("tibble")

#############################################################
#set working directory
setwd("C:/Users/catar/Documents/Dundee PhD/Thesis/Chapter 3 Analysis")
getwd()

#Import the datasets
pgpt <- read_tsv("sample_by_pgpt_family_presence_relative_isolate_abundance.tsv")
KO <- read_tsv("sample_by_ko_presence_relative_isolate_abundance.tsv")

#Inspect the files
View(pgpt)
View(KO)

#Convert to matrix
pgpt <- pgpt %>%
  tibble::column_to_rownames("sample_id") 
pgpt <- as.matrix(pgpt)
pgpt <- t(pgpt)
class(pgpt)

KO <- KO %>%
  tibble::column_to_rownames("sample_id") 
KO <- as.matrix(KO)
KO <- t(KO)
class(KO)

#Create taxonomy matrix
PGPTtaxmat = matrix(rownames(pgpt))
rownames(PGPTtaxmat) <- rownames(pgpt)
colnames(PGPTtaxmat) <- c("Gene_ID")
View(PGPTtaxmat)
class(PGPTtaxmat)

KOtaxmat = matrix(rownames(KO))
rownames(KOtaxmat) <- rownames(KO)
colnames(KOtaxmat) <- c("KO_Term")
View(KOtaxmat)
class(KOtaxmat)

#Mapping file
library(readxl)
Sample_Map <- read_excel("Sample_Map.xlsx")
View(Sample_Map)

#Set rownames
Sample_Map <- Sample_Map %>% 
  tibble::column_to_rownames("Sample_ID") 
View(Sample_Map)

#Create phyloseq objects
OTU = otu_table(pgpt, taxa_are_rows = TRUE)
TAX = tax_table(PGPTtaxmat)
samples = sample_data(Sample_Map)
pgptPhylo <- phyloseq(OTU, TAX, samples)
pgptPhylo

OTU = otu_table(KO, taxa_are_rows = TRUE)
TAX = tax_table(KOtaxmat)
samples = sample_data(Sample_Map)
koPhylo <- phyloseq(OTU, TAX, samples)
koPhylo

#######################
####Start with pgpt####
#######################

#subset for KAS1
pgpt_KAS1 <- subset_samples(pgptPhylo, SynCom_ID == "KAS1")
pgpt_KAS1
sample_names(pgpt_KAS1)

#Any genes with no abundance?
pgpt_KAS1_prune <- prune_taxa(taxa_sums(pgpt_KAS1) > 0, pgpt_KAS1)
pgpt_KAS1_prune

#Order factors
sample_data(pgpt_KAS1)$Genotype <- factor(sample_data(pgpt_KAS1)$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
sample_data(pgpt_KAS1)$Plant <- factor(sample_data(pgpt_KAS1)$Plant, levels=c("Unplanted", "Planted"))
sample_data(pgpt_KAS1)$Date_Inoc <- as.factor(sample_data(pgpt_KAS1)$Date_Inoc)

#Set colours
SynCom_col <- c("#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

####Whole set####

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
pgpt_KAS1.ord <- ordinate(pgpt_KAS1, "NMDS", "bray")
#plot
p = plot_ordination(pgpt_KAS1, pgpt_KAS1.ord, type="samples", color="Genotype", shape="Syncom_Status") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#Constrained for SynCom status
pgpt_KAS1.cap <- ordinate(pgpt_KAS1, "CAP", "bray", ~ Syncom_Status)
#Plot
a = plot_ordination(pgpt_KAS1, pgpt_KAS1.cap, type="samples", color="Genotype", shape="Syncom_Status") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Constrained for genotype
pgpt_KAS1.cap <- ordinate(pgpt_KAS1, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(pgpt_KAS1, pgpt_KAS1.cap, type="samples", color="Genotype", shape="Syncom_Status") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Stats - Permanova
BC_bacteria  <- phyloseq::distance(pgpt_KAS1, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Plant * Syncom_Status * Date_Inoc, data= as.data.frame(as.matrix(sample_data(pgpt_KAS1))), permutations = 5000)
Stat

####Live Only####
#subset for live
pgpt_KAS1_live <- subset_samples(pgpt_KAS1, Syncom_Status == "Live")
pgpt_KAS1_live

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
pgpt_KAS1_live.ord <- ordinate(pgpt_KAS1_live, "NMDS", "bray")
#plot
p = plot_ordination(pgpt_KAS1_live, pgpt_KAS1_live.ord, type="samples", color="Genotype", shape="Date_Inoc") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#Constrained for run
pgpt_KAS1_live.cap <- ordinate(pgpt_KAS1_live, "CAP", "bray", ~ Date_Inoc)
#Plot
a = plot_ordination(pgpt_KAS1_live, pgpt_KAS1_live.cap, type="samples", color="Genotype", shape="Date_Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Constrained for genotype
pgpt_KAS1_live.cap <- ordinate(pgpt_KAS1_live, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(pgpt_KAS1_live, pgpt_KAS1_live.cap, type="samples", color="Genotype", shape="Date_Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Stats - Permanova
BC_bacteria  <- phyloseq::distance(pgpt_KAS1_live, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Plant * Date_Inoc, data= as.data.frame(as.matrix(sample_data(pgpt_KAS1_live))), permutations = 5000)
Stat

#Deseq
#extract count data 
pgpt_KAS1_live_ex <- otu_table(pgpt_KAS1_live)
countData = as.data.frame(pgpt_KAS1_live_ex)
countData <- countData*1000000
countData <- round(countData)
colnames(pgpt_KAS1_live_ex)
pgpt_KAS1_live_ex <- pgpt_KAS1_live_ex*1000000
pgpt_KAS1_live_ex <- round(pgpt_KAS1_live_ex)
pgpt_KAS1_live_ex <- as.integer(pgpt_KAS1_live_ex)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(pgpt_KAS1_live)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
pgpt_KAS1_live_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Plant)

#execute the differential count analysis with the function DESeq 
pgpt_KAS1_live_test <- DESeq(pgpt_KAS1_live_cds, fitType="local", betaPrior = FALSE) 

#define the Genera differentially enriched in Barke vs Bulk
pgpt_KAS1_live_test_plant <- results(pgpt_KAS1_live_test, contrast = c("Plant", "Unplanted", "Planted")) 

#extract  ASVs whose adjusted p.value in a given comparison is below 0.05
pgpt_KAS1_live_test_plant_sig <- pgpt_KAS1_live_test_plant[(rownames(pgpt_KAS1_live_test_plant)[which(pgpt_KAS1_live_test_plant$padj <0.05)]), ]

#inspect the generated file
summary(pgpt_KAS1_live_test_plant_sig)

up <- pgpt_KAS1_live_test_plant_sig[which(pgpt_KAS1_live_test_plant_sig$log2FoldChange > 0.585 & pgpt_KAS1_live_test_plant_sig$padj < 0.05),]
down <- pgpt_KAS1_live_test_plant_sig[which(pgpt_KAS1_live_test_plant_sig$log2FoldChange < -0.585 & pgpt_KAS1_live_test_plant_sig$padj < 0.05),]

#write.csv(up, file="PGPT_Unplanted.csv")
#write.csv(down, file="PGPT_Planted.csv")

####HI Only####
#subset for HI
pgpt_KAS1_HI <- subset_samples(pgpt_KAS1, Syncom_Status == "HK")
pgpt_KAS1_HI

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
pgpt_KAS1_HI.ord <- ordinate(pgpt_KAS1_HI, "NMDS", "bray")
#plot
p = plot_ordination(pgpt_KAS1_HI, pgpt_KAS1_HI.ord, type="samples", color="Genotype", shape="Date_Inoc") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#Constrained for run
pgpt_KAS1_HI.cap <- ordinate(pgpt_KAS1_HI, "CAP", "bray", ~ Date_Inoc)
#Plot
a = plot_ordination(pgpt_KAS1_HI, pgpt_KAS1_HI.cap, type="samples", color="Genotype", shape="Date_Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Constrained for genotype
pgpt_KAS1_HI.cap <- ordinate(pgpt_KAS1_HI, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(pgpt_KAS1_HI, pgpt_KAS1_HI.cap, type="samples", color="Genotype", shape="Date_Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Stats - Permanova
BC_bacteria  <- phyloseq::distance(pgpt_KAS1_HI, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Plant * Date_Inoc, data= as.data.frame(as.matrix(sample_data(pgpt_KAS1_HI))), permutations = 5000)
Stat

####Planted only####
#subset for plant
pgpt_KAS1_plant <- subset_samples(pgpt_KAS1, Plant == "Planted")
pgpt_KAS1_plant

#set the colours: planted
SynCom_col_p <- c("#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
pgpt_KAS1_plant.ord <- ordinate(pgpt_KAS1_plant, "NMDS", "bray")
#plot
p = plot_ordination(pgpt_KAS1_plant, pgpt_KAS1_plant.ord, type="samples", color="Genotype", shape="Syncom_Status") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col_p)
p

#Constrained for SynCom status
pgpt_KAS1_plant.cap <- ordinate(pgpt_KAS1_plant, "CAP", "bray", ~ Syncom_Status)
#Plot
a = plot_ordination(pgpt_KAS1_plant, pgpt_KAS1_plant.cap, type="samples", color="Genotype", shape="Syncom_Status") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_p)
a

#Constrained for genotype
pgpt_KAS1_plant.cap <- ordinate(pgpt_KAS1_plant, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(pgpt_KAS1_plant, pgpt_KAS1_plant.cap, type="samples", color="Genotype", shape="Syncom_Status") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_p)
a

#Stats - Permanova
BC_bacteria  <- phyloseq::distance(pgpt_KAS1_plant, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Genotype * Syncom_Status * Date_Inoc, data= as.data.frame(as.matrix(sample_data(pgpt_KAS1_plant))), permutations = 5000)
Stat

####Live Only####
#subset for live
pgpt_KAS1_plant_live <- subset_samples(pgpt_KAS1_plant, Syncom_Status == "Live")
pgpt_KAS1_plant_live

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
pgpt_KAS1_plant_live.ord <- ordinate(pgpt_KAS1_plant_live, "NMDS", "bray")
#plot
p = plot_ordination(pgpt_KAS1_plant_live, pgpt_KAS1_plant_live.ord, type="samples", color="Genotype", shape="Date_Inoc") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col_p)
p

#Constrained for run
pgpt_KAS1_plant_live.cap <- ordinate(pgpt_KAS1_plant_live, "CAP", "bray", ~ Date_Inoc)
#Plot
a = plot_ordination(pgpt_KAS1_plant_live, pgpt_KAS1_plant_live.cap, type="samples", color="Genotype", shape="Date_Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_p)
a

#Constrained for genotype
pgpt_KAS1_plant_live.cap <- ordinate(pgpt_KAS1_plant_live, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(pgpt_KAS1_plant_live, pgpt_KAS1_plant_live.cap, type="samples", color="Genotype", shape="Date_Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_p)
a

#Stats - Permanova
BC_bacteria  <- phyloseq::distance(pgpt_KAS1_plant_live, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Genotype * Date_Inoc, data= as.data.frame(as.matrix(sample_data(pgpt_KAS1_plant_live))), permutations = 5000)
Stat

####HI Only####
#subset for HI
pgpt_KAS1_plant_HI <- subset_samples(pgpt_KAS1_plant, Syncom_Status == "HK")
pgpt_KAS1_plant_HI

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
pgpt_KAS1_plant_HI.ord <- ordinate(pgpt_KAS1_plant_HI, "NMDS", "bray")
#plot
p = plot_ordination(pgpt_KAS1_plant_HI, pgpt_KAS1_plant_HI.ord, type="samples", color="Genotype", shape="Date_Inoc") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col_p)
p

#Constrained for run
pgpt_KAS1_plant_HI.cap <- ordinate(pgpt_KAS1_plant_HI, "CAP", "bray", ~ Date_Inoc)
#Plot
a = plot_ordination(pgpt_KAS1_plant_HI, pgpt_KAS1_plant_HI.cap, type="samples", color="Genotype", shape="Date_Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_p)
a

#Constrained for genotype
pgpt_KAS1_plant_HI.cap <- ordinate(pgpt_KAS1_plant_HI, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(pgpt_KAS1_plant_HI, pgpt_KAS1_plant_HI.cap, type="samples", color="Genotype", shape="Date_Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_p)
a

#Stats - Permanova
BC_bacteria  <- phyloseq::distance(pgpt_KAS1_plant_HI, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Genotype * Date_Inoc, data= as.data.frame(as.matrix(sample_data(pgpt_KAS1_plant_HI))), permutations = 5000)
Stat

#####################
####Move on to KO####
#####################

#subset for KAS1
ko_KAS1 <- subset_samples(koPhylo, SynCom_ID == "KAS1")
ko_KAS1

#Order factors
sample_data(ko_KAS1)$Genotype <- factor(sample_data(ko_KAS1)$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
sample_data(ko_KAS1)$Plant <- factor(sample_data(ko_KAS1)$Plant, levels=c("Unplanted", "Planted"))
sample_data(ko_KAS1)$Date_Inoc <- as.factor(sample_data(ko_KAS1)$Date_Inoc)

#Set colours
SynCom_col <- c("#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

####Whole set####

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
ko_KAS1.ord <- ordinate(ko_KAS1, "NMDS", "bray")
#plot
p = plot_ordination(ko_KAS1, ko_KAS1.ord, type="samples", color="Genotype", shape="Syncom_Status") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

####Thesis####
#Constrained for SynCom status
ko_KAS1.cap <- ordinate(ko_KAS1, "CAP", "bray", ~ Syncom_Status)
#Plot
a = plot_ordination(ko_KAS1, ko_KAS1.cap, type="samples", color="Genotype", shape="Syncom_Status") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)+
  scale_shape_discrete(name = "SynCom Status", labels = c("HI", "Live")) + 
  theme(legend.position="bottom") 
a

#Constrained for genotype
ko_KAS1.cap <- ordinate(ko_KAS1, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(ko_KAS1, ko_KAS1.cap, type="samples", color="Genotype", shape="Syncom_Status") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Stats - Permanova
BC_bacteria  <- phyloseq::distance(ko_KAS1, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Plant * Syncom_Status * Date_Inoc, data= as.data.frame(as.matrix(sample_data(ko_KAS1))), permutations = 5000)
Stat

####Live Only####
#subset for live
ko_KAS1_live <- subset_samples(ko_KAS1, Syncom_Status == "Live")
ko_KAS1_live

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
ko_KAS1_live.ord <- ordinate(ko_KAS1_live, "NMDS", "bray")
#plot
p = plot_ordination(ko_KAS1_live, ko_KAS1_live.ord, type="samples", color="Genotype", shape="Date_Inoc") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#Thesis
#Constrained for run
ko_KAS1_live.cap <- ordinate(ko_KAS1_live, "CAP", "bray", ~ Date_Inoc)
#Plot
a = plot_ordination(ko_KAS1_live, ko_KAS1_live.cap, type="samples", color="Genotype", shape="Date_Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)+
  scale_shape_discrete(name = "Replicate", labels = c("1", "2")) + 
  theme(legend.position="bottom") 
a

#Constrained for genotype
ko_KAS1_live.cap <- ordinate(ko_KAS1_live, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(ko_KAS1_live, ko_KAS1_live.cap, type="samples", color="Genotype", shape="Date_Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Stats - Permanova
BC_bacteria  <- phyloseq::distance(ko_KAS1_live, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Plant * Date_Inoc, data= as.data.frame(as.matrix(sample_data(ko_KAS1_live))), permutations = 5000)
Stat

#Deseq
#extract count data 
ko_KAS1_live_ex <- otu_table(ko_KAS1_live)
countData = as.data.frame(ko_KAS1_live_ex)
countData <- countData*1000000
countData <- round(countData)
colnames(ko_KAS1_live_ex)
ko_KAS1_live_ex <- ko_KAS1_live_ex*1000000
ko_KAS1_live_ex <- round(ko_KAS1_live_ex)
ko_KAS1_live_ex <- as.integer(ko_KAS1_live_ex)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(ko_KAS1_live)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
ko_KAS1_live_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Plant)

#execute the differential count analysis with the function DESeq 
ko_KAS1_live_test <- DESeq(ko_KAS1_live_cds, fitType="local", betaPrior = FALSE) 

#define the Genera differentially enriched in Barke vs Bulk
ko_KAS1_live_test_plant <- results(ko_KAS1_live_test, contrast = c("Plant", "Unplanted", "Planted")) 

#extract  ASVs whose adjusted p.value in a given comparison is below 0.05
ko_KAS1_live_test_plant_sig <- ko_KAS1_live_test_plant[(rownames(ko_KAS1_live_test_plant)[which(ko_KAS1_live_test_plant$padj <0.05)]), ]

#inspect the generated file
summary(ko_KAS1_live_test_plant_sig)

up <- ko_KAS1_live_test_plant_sig[which(ko_KAS1_live_test_plant_sig$log2FoldChange > 0.585 & ko_KAS1_live_test_plant_sig$padj < 0.05),]
down <- ko_KAS1_live_test_plant_sig[which(ko_KAS1_live_test_plant_sig$log2FoldChange < -0.585 & ko_KAS1_live_test_plant_sig$padj < 0.05),]

#write.csv(up, file="ko_Unplanted.csv")
#write.csv(down, file="ko_Planted.csv")

####HI Only####
#subset for HI
ko_KAS1_HI <- subset_samples(ko_KAS1, Syncom_Status == "HK")
ko_KAS1_HI

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
ko_KAS1_HI.ord <- ordinate(ko_KAS1_HI, "NMDS", "bray")
#plot
p = plot_ordination(ko_KAS1_HI, ko_KAS1_HI.ord, type="samples", color="Genotype", shape="Date_Inoc") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

####Thesis####
#Constrained for run
ko_KAS1_HI.cap <- ordinate(ko_KAS1_HI, "CAP", "bray", ~ Date_Inoc)
#Plot
a = plot_ordination(ko_KAS1_HI, ko_KAS1_HI.cap, type="samples", color="Genotype", shape="Date_Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)+
  scale_shape_discrete(name = "Replicate", labels = c("1", "2")) + 
  theme(legend.position="bottom") 
a

#Constrained for genotype
ko_KAS1_HI.cap <- ordinate(ko_KAS1_HI, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(ko_KAS1_HI, ko_KAS1_HI.cap, type="samples", color="Genotype", shape="Date_Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Stats - Permanova
BC_bacteria  <- phyloseq::distance(ko_KAS1_HI, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Plant * Date_Inoc, data= as.data.frame(as.matrix(sample_data(ko_KAS1_HI))), permutations = 5000)
Stat

####Planted only####
#subset for plant
ko_KAS1_plant <- subset_samples(ko_KAS1, Plant == "Planted")
ko_KAS1_plant

#set the colours: planted
SynCom_col_p <- c("#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
ko_KAS1_plant.ord <- ordinate(ko_KAS1_plant, "NMDS", "bray")
#plot
p = plot_ordination(ko_KAS1_plant, ko_KAS1_plant.ord, type="samples", color="Genotype", shape="Syncom_Status") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col_p)
p

####Thesis####
#Constrained for SynCom status
ko_KAS1_plant.cap <- ordinate(ko_KAS1_plant, "CAP", "bray", ~ Syncom_Status)
#Plot
a = plot_ordination(ko_KAS1_plant, ko_KAS1_plant.cap, type="samples", color="Genotype", shape="Syncom_Status") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_p)+
  scale_shape_discrete(name = "SynCom Status", labels = c("HI", "Live")) + 
  theme(legend.position="bottom") 
a

#Constrained for genotype
ko_KAS1_plant.cap <- ordinate(ko_KAS1_plant, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(ko_KAS1_plant, ko_KAS1_plant.cap, type="samples", color="Genotype", shape="Syncom_Status") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_p)
a

#Stats - Permanova
BC_bacteria  <- phyloseq::distance(ko_KAS1_plant, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Genotype * Syncom_Status * Date_Inoc, data= as.data.frame(as.matrix(sample_data(ko_KAS1_plant))), permutations = 5000)
Stat

####Live Only####
#subset for live
ko_KAS1_plant_live <- subset_samples(ko_KAS1_plant, Syncom_Status == "Live")
ko_KAS1_plant_live

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
ko_KAS1_plant_live.ord <- ordinate(ko_KAS1_plant_live, "NMDS", "bray")
#plot
p = plot_ordination(ko_KAS1_plant_live, ko_KAS1_plant_live.ord, type="samples", color="Genotype", shape="Date_Inoc") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col_p)
p

#Constrained for run
ko_KAS1_plant_live.cap <- ordinate(ko_KAS1_plant_live, "CAP", "bray", ~ Date_Inoc)
#Plot
a = plot_ordination(ko_KAS1_plant_live, ko_KAS1_plant_live.cap, type="samples", color="Genotype", shape="Date_Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_p)
a

#Constrained for genotype
ko_KAS1_plant_live.cap <- ordinate(ko_KAS1_plant_live, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(ko_KAS1_plant_live, ko_KAS1_plant_live.cap, type="samples", color="Genotype", shape="Date_Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_p)
a

#Stats - Permanova
BC_bacteria  <- phyloseq::distance(ko_KAS1_plant_live, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Genotype * Date_Inoc, data= as.data.frame(as.matrix(sample_data(ko_KAS1_plant_live))), permutations = 5000)
Stat

####HI Only####
#subset for HI
ko_KAS1_plant_HI <- subset_samples(ko_KAS1_plant, Syncom_Status == "HK")
ko_KAS1_plant_HI

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
ko_KAS1_plant_HI.ord <- ordinate(ko_KAS1_plant_HI, "NMDS", "bray")
#plot
p = plot_ordination(ko_KAS1_plant_HI, ko_KAS1_plant_HI.ord, type="samples", color="Genotype", shape="Date_Inoc") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col_p)
p

#Constrained for run
ko_KAS1_plant_HI.cap <- ordinate(ko_KAS1_plant_HI, "CAP", "bray", ~ Date_Inoc)
#Plot
a = plot_ordination(ko_KAS1_plant_HI, ko_KAS1_plant_HI.cap, type="samples", color="Genotype", shape="Date_Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_p)
a

#Constrained for genotype
ko_KAS1_plant_HI.cap <- ordinate(ko_KAS1_plant_HI, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(ko_KAS1_plant_HI, ko_KAS1_plant_HI.cap, type="samples", color="Genotype", shape="Date_Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_p)
a

#Stats - Permanova
BC_bacteria  <- phyloseq::distance(ko_KAS1_plant_HI, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Genotype * Date_Inoc, data= as.data.frame(as.matrix(sample_data(ko_KAS1_plant_HI))), permutations = 5000)
Stat
