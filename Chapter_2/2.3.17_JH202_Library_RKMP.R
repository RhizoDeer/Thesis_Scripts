#####################################################################################
#PREPROCESING
####################################################################################
#############################################################
#
# Ref to the ARTICLE 
# 
#  Code to compute calculations SynComs RK & KA
#  Revision 05/25
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
#required packages 
library("phyloseq")
library("vegan")
library ("ggplot2")
library ("DESeq2")
library("viridis")
library("gridExtra")

#############################################################
#set working directory: DB cpu
#setwd("/cluster/db/R_shared/NSQ04_v0525/")
#set working directory: Katie home
setwd("C:/Users/catar/Documents/Dundee PhD/SynCom/JH202")
#set working dirctory: Katie work
#setwd("~/PhD Files/SynCom/JH202")

#import the datasets
NSQ04 <-readRDS("NSQ04_dada2_silva_138.1.rds")

#inspect the files
NSQ04

#total number of reads
sum(sample_sums(NSQ04))

#and its three "constituents"

#ASV counts
otu_table(NSQ04)

#Taxonomy information
tax_table(NSQ04)

#mapping files
sample_data(NSQ04)

#stats overall library
sum(sample_sums(NSQ04))

#proceed with independent filtering in both objects prior merging
##################################################################
#Pre-processing: remove Chloroplast and Mitochondria but retain NA
#rationale: 16S rRNA primers may amplify plant-derived sequences 
#those may interfere with data analysis (we are after bacteria, not plants)
#################################################################

#NSQ04
NSQ04_no_chlor <-subset_taxa(NSQ04, (Order!="Chloroplast") | is.na(Order))
NSQ04_no_chlor

NSQ04_no_plants <-subset_taxa(NSQ04_no_chlor, (Family!="Mitochondria") | is.na(Family))
NSQ04_no_plants

##################################################################
#Prune putative contaminant ASVs
#Rationale: bacteria are everywhere and 16S rRNA amplification is very efficient
#We use a list of bacteria likely representing the contamination we may face in our lab as baseline
#More info here: https://www.frontiersin.org/articles/10.3389/fmicb.2018.01650/full 
##################################################################

#Import the list of contaminat ASVs from JH06 library
Contaminant_ASVs <- read.delim("JH06_contaminant_ASVs_ids.txt", header = FALSE)

#identify the proportion of putative contaminants in the merged object
NSQ04_contaminants <- intersect(taxa_names(NSQ04_no_plants),Contaminant_ASVs)
NSQ04_contaminants

#########################################################################
# Remove ASVs assigned to NA at phylum level
# Rationale: if a sequence cannot be assigned at a very high level such as Phylum is not of use us
#########################################################################

NSQ04_no_plants_1 <- subset_taxa(NSQ04_no_plants, Phylum!= "NA")
NSQ04_no_plants_1

#########################################################################
#Subsest for RKMP
#remove editing experiment
#remove ASV with 0 counts
#########################################################################

NSQ04_RKMP <- subset_samples(NSQ04_no_plants_1, SynCom.ID != "KAS1")
NSQ04_RKMP_integer <- prune_taxa(taxa_sums(NSQ04_RKMP) > 0, NSQ04_RKMP)
NSQ04_RKMP_integer

#########################################################################
# inspect reads distribution across objects
# Rationale: if we want to compare like-with-like we should know sample distribution
#########################################################################

#summary
sum(sample_sums(NSQ04_RKMP))
min(sample_sums(NSQ04_RKMP))
max(sample_sums(NSQ04_RKMP))
mean(sample_sums(NSQ04_RKMP))

#graphical outputs
sort(sample_sums(NSQ04_RKMP))
hist(sample_sums(NSQ04_RKMP), main = paste("NSQ04 reads distribution; RKMP"), xlab = paste("reads"), ylab = paste ("number of samples"))

#inspect data distribution of various experiments: Heat Killed
NSQ04_RKMP_HK <- subset_samples(NSQ04_RKMP, Syncom.Status == "HK")
sort(sample_sums(NSQ04_RKMP_HK))
hist(sample_sums(NSQ04_RKMP_HK), main = paste("NSQ04 reads distribution; RKMP_HK"), xlab = paste("reads"), ylab = paste ("number of samples"))

#inspect data distribution of various experiments: Live
NSQ04_RKMP_Live <- subset_samples(NSQ04_RKMP, Syncom.Status == "Live")
sort(sample_sums(NSQ04_RKMP_Live))
hist(sample_sums(NSQ04_RKMP_Live), main = paste("NSQ04 reads distribution; RKMP_Live"), xlab = paste("reads"), ylab = paste ("number of samples"))

#inspect data distribution of various experiments: inoculum
NSQ04_RKMP_inoculum <- subset_samples(NSQ04_RKMP, Genotype == "Inoculum")
sort(sample_sums(NSQ04_RKMP_inoculum))
hist(sample_sums(NSQ04_RKMP_inoculum), main = paste("NSQ04 reads distribution; RKMP_Inoculum"), xlab = paste("reads"), ylab = paste ("number of samples"))

#inspect data distribution of various experiments: PBS controls
NSQ04_RKMP_PBS <- subset_samples(NSQ04_RKMP, SynCom.ID == "PBS")
sort(sample_sums(NSQ04_RKMP_PBS))
hist(sample_sums(NSQ04_RKMP_PBS), main = paste("NSQ04 reads distribution; RKMP_PBS"), xlab = paste("reads"), ylab = paste ("number of samples"))

#########################################################################
# work with inoculum samples and set criteria to work with unevenly sequenced samples
# a) keep ASV with more than 20 reads in at least 10% of the samples
# b) rarefy at even sequencing depth
# c) agglomerate at genus level
#########################################################################

#set an arbitrary threshold of 20 reads in 10% of samples
NSQ04_RKMP_inoculum_threshold = filter_taxa(NSQ04_RKMP_inoculum, function(x) sum(x > 20) > (0.1 *length(x)), TRUE)
NSQ04_RKMP_inoculum_threshold

#inspect read distribution upon thresholding
sort(sample_sums(NSQ04_RKMP_inoculum_threshold))

#agglomerate at genus level
NSQ04_RKMP_inoc_genus <- tax_glom(NSQ04_RKMP_inoculum_threshold, taxrank="Genus")

#inspect the files
NSQ04_RKMP_inoc_genus
sort(sample_sums(NSQ04_RKMP_inoc_genus))

#rarefy at 34K
NSQ04_RKMP_inoc_threshold_34K <- rarefy_even_depth(NSQ04_RKMP_inoc_genus, 34000)

#data visualisation
p = plot_bar(NSQ04_RKMP_inoc_threshold_34K, "Genus", facet_grid=Syncom.Status ~ SynCom.ID, fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#save the file for the reproducibility of the code
#saveRDS(NSQ04_RKMP_inoc_threshold_34K, file = "NSQ04_RKMP_inoc_calculation.rds")
NSQ04_RKMP_inoc_calculation <-readRDS("NSQ04_RKMP_inoc_calculation.rds")
NSQ04_RKMP_inoc_calculation

#########################################################################
# work with live samples and set criteria to work with unevenly sequenced samples
# a) remove samples with less than 1,000 reads
# b) remove ASV with less than 15 reads in at least 10% of the samples
# c) rarefy at even sequencing depth
# d) agglomerate at genus level
#########################################################################

#remove poorly sequenced samples
NSQ04_RKMP_Live_19K <- prune_samples(sample_sums(NSQ04_RKMP_Live) > 19000, NSQ04_RKMP_Live)

#set an arbitrary threshold of 20 reads in 10% of samples
NSQ04_RKMP_Live_19K_threshold = filter_taxa(NSQ04_RKMP_Live_19K, function(x) sum(x > 20) > (0.1 *length(x)), TRUE)
NSQ04_RKMP_Live_19K_threshold

#inspect read distribution upon thresholding
sort(sample_sums(NSQ04_RKMP_Live_19K_threshold))

#agglomerate at genus level
NSQ04_RKMP_Live_genus <- tax_glom(NSQ04_RKMP_Live_19K_threshold, taxrank="Genus")

#inspect the files
NSQ04_RKMP_Live_genus
sort(sample_sums(NSQ04_RKMP_Live_genus))

#rarefy at 19K
NSQ04_RKMP_Live_threshold_19K <- rarefy_even_depth(NSQ04_RKMP_Live_genus, 19000)

#data visualisation
p = plot_bar(NSQ04_RKMP_Live_threshold_19K, "Genus", facet_grid=SynCom.ID ~ Genotype, fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#save the file for the reproducibility of the code
#saveRDS(NSQ04_RKMP_Live_threshold_19K, file = "NSQ04_RKMP_Live_calculation.rds")
NSQ04_RKMP_calculation <-readRDS("NSQ04_RKMP_Live_calculation.rds")
NSQ04_RKMP_calculation

##############################################################
#All live
##############################################################

#re-order the factors
sample_data(NSQ04_RKMP_calculation)$Genotype <- factor(sample_data(NSQ04_RKMP_calculation)$Genotype, levels=c("Inoculum", "Bulk", "Barke", "124-17", "124-52", "Morex"))
sample_data(NSQ04_RKMP_calculation)$SynCom.ID <- factor(sample_data(NSQ04_RKMP_calculation)$SynCom.ID, levels=c("RKMP.Old", "RKMP.New"))
sample_data(NSQ04_RKMP_calculation)$Plant <- factor(sample_data(NSQ04_RKMP_calculation)$Plant, levels=c("Unplanted", "Planted"))
sample_data(NSQ04_RKMP_calculation)$Date.Inoc <- as.factor(sample_data(NSQ04_RKMP_calculation)$Date.Inoc)

#set the colours: with inoculum
SynCom_col <- c("#999999", "#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#Stacked bar plot
p = plot_bar(NSQ04_RKMP_calculation, "Genus", facet_grid=SynCom.ID ~ Genotype, fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#alpha diversity calculation
#https://www.nature.com/articles/s41598-024-77864-y
p = plot_richness(NSQ04_RKMP_calculation, x="Genotype", color="Genotype", measures=c("Observed", "Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#remove Observed, this is just a confirmation
p = plot_richness(NSQ04_RKMP_calculation, x="Genotype", color="Genotype", measures=c("Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#remove inoculum for plotting and statistical analysis
NSQ04_RKMP_calculation_ni <- subset_samples(NSQ04_RKMP_calculation, Genotype != "Inoculum")

#set the colours: no inoculum
SynCom_col_ni <- c("#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#without inoculum
p = plot_richness(NSQ04_RKMP_calculation_ni, x="Genotype", color="Genotype", measures=c("Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col_ni)
p

#ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
NSQ04_RKMP.ord <- ordinate(NSQ04_RKMP_calculation, "NMDS", "bray")

#ggplots function to increase effectivness of the visualisation
p = plot_ordination(NSQ04_RKMP_calculation, NSQ04_RKMP.ord, type="samples", color="Genotype", shape="SynCom.ID") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#alternative plotting were we "constrain for" for genotype effect
NSQ04_RKMP_G.cap <- ordinate(NSQ04_RKMP_calculation, "CAP", "bray", ~ Genotype)

#ggplots function to increase effectivness of the visualisation
a = plot_ordination(NSQ04_RKMP_calculation, NSQ04_RKMP_G.cap, type="samples", color="Genotype", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#alternative plotting were we "constrain for" for SynCom effect
NSQ04_RKMP_S.cap <- ordinate(NSQ04_RKMP_calculation, "CAP", "bray", ~ SynCom.ID)

#ggplots function to increase effectivness of the visualisation
a = plot_ordination(NSQ04_RKMP_calculation, NSQ04_RKMP_S.cap, type="samples", color="Genotype", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

##############################################################
#Run 1 Only
##############################################################

#re-order the factors
sample_data(NSQ04_RKMP_calculation)$Genotype <- factor(sample_data(NSQ04_RKMP_calculation)$Genotype, levels=c("Inoculum", "Bulk", "Barke", "124-17", "124-52", "Morex"))
sample_data(NSQ04_RKMP_calculation)$SynCom.ID <- factor(sample_data(NSQ04_RKMP_calculation)$SynCom.ID, levels=c("RKMP.Old", "RKMP.New"))
sample_data(NSQ04_RKMP_calculation)$Plant <- factor(sample_data(NSQ04_RKMP_calculation)$Plant, levels=c("Unplanted", "Planted"))
sample_data(NSQ04_RKMP_calculation)$Date.Inoc <- as.factor(sample_data(NSQ04_RKMP_calculation)$Date.Inoc)

#set the colours: with inoculum
SynCom_col <- c("#999999", "#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")
#set the colours: no inoculum
SynCom_col_ni <- c("#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")
#set the colours: plant
Plant_col <- c("#E1BE6A", "#40B0A6")

#Subset for run 1
NSQ04_RKMP_R1 <- subset_samples(NSQ04_RKMP_calculation, Date.Inoc == "241003")
sort(sample_sums(NSQ04_RKMP_R1))

#Stacked bar plot
p = plot_bar(NSQ04_RKMP_R1, "Genus", facet_grid=SynCom.ID ~ Genotype, fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#Without inoculum
NSQ04_RKMP_R1_ni <- subset_samples(NSQ04_RKMP_R1, Genotype != "Inoculum")
sort(sample_sums(NSQ04_RKMP_R1_ni))

#Stacked bar plot
p = plot_bar(NSQ04_RKMP_R1_ni, "Genus", facet_grid=SynCom.ID ~ Plant, fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#alpha diversity calculation
#https://www.nature.com/articles/s41598-024-77864-y
p = plot_richness(NSQ04_RKMP_R1, x="SynCom.ID", color="Genotype", measures=c("Observed", "Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#remove Observed, this is just a confirmation
p = plot_richness(NSQ04_RKMP_R1, x="SynCom.ID", color="Genotype", measures=c("Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#without inoculum
p = plot_richness(NSQ04_RKMP_R1_ni, x="SynCom.ID", color="Plant", measures=c("Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = Plant_col)
p

####PAPER####
#without inoculum but coloured by genotype and with observed
p = plot_richness(NSQ04_RKMP_R1_ni, x="SynCom.ID", color="Genotype", measures=c("Observed","Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col_ni)+
  scale_x_discrete(labels = c("By OD", "By CFU"))+
  xlab("")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "bottom", axis.title.x = element_blank())
p

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
NSQ04_RKMP_R1.ord <- ordinate(NSQ04_RKMP_R1, "NMDS", "bray")
#plot
p = plot_ordination(NSQ04_RKMP_R1, NSQ04_RKMP_R1.ord, type="samples", color="Genotype", shape="SynCom.ID") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#unconstrained ordination - no inoculum
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
NSQ04_RKMP_R1_ni.ord <- ordinate(NSQ04_RKMP_R1_ni, "NMDS", "bray")
#plot
p = plot_ordination(NSQ04_RKMP_R1_ni, NSQ04_RKMP_R1_ni.ord, type="samples", color="Plant", shape="SynCom.ID") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = Plant_col)
p

#Constrained for SynCom effect
NSQ04_RKMP_R1_S.cap <- ordinate(NSQ04_RKMP_R1, "CAP", "bray", ~ SynCom.ID)
#Plot
a = plot_ordination(NSQ04_RKMP_R1, NSQ04_RKMP_R1_S.cap, type="samples", color="Genotype", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Constrained for SynCom effect - no inoculum
NSQ04_RKMP_R1_ni_S.cap <- ordinate(NSQ04_RKMP_R1_ni, "CAP", "bray", ~ SynCom.ID)
#Plot
a = plot_ordination(NSQ04_RKMP_R1_ni, NSQ04_RKMP_R1_ni_S.cap, type="samples", color="Plant", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = Plant_col)
a

#Constrained for 'Genotype' effect
NSQ04_RKMP_R1_G.cap <- ordinate(NSQ04_RKMP_R1, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(NSQ04_RKMP_R1, NSQ04_RKMP_R1_G.cap, type="samples", color="Genotype", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Constrained for plant effect - no inoculum
NSQ04_RKMP_R1_ni_P.cap <- ordinate(NSQ04_RKMP_R1_ni, "CAP", "bray", ~ Plant)
#Plot
a = plot_ordination(NSQ04_RKMP_R1_ni, NSQ04_RKMP_R1_ni_P.cap, type="samples", color="Plant", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = Plant_col)
a

#Stats - Permanova
BC_bacteria  <- phyloseq::distance(NSQ04_RKMP_R1_ni, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Plant * SynCom.ID, data= as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R1_ni))), permutations = 5000)
Stat

#Deseq
#extract count data 
NSQ04_RKMP_genotype_R1_ni <- otu_table(NSQ04_RKMP_R1_ni)
countData = as.data.frame(NSQ04_RKMP_genotype_R1_ni)
colnames(NSQ04_RKMP_genotype_R1_ni)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R1_ni)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_RKMP_R1_ni_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)

#execute the differential count analysis with the function DESeq 
NSQ04_RKMP_R1_ni_test <- DESeq(NSQ04_RKMP_R1_ni_cds, fitType="local", betaPrior = FALSE) 

#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_R1_ni_Ba <- results(NSQ04_RKMP_R1_ni_test, contrast = c("Genotype", "Bulk", "Barke")) 

#extract  ASVs whose adjusted p.value in a given comparison is below 0.05 [if any]
rhizo_genotype_RKMP_R1_ni_Ba_FDR005 <- rhizo_genotype_RKMP_R1_ni_Ba[(rownames(rhizo_genotype_RKMP_R1_ni_Ba)[which(rhizo_genotype_RKMP_R1_ni_Ba$padj <0.05)]), ]

#inspect the generated file
rhizo_genotype_RKMP_R1_ni_Ba_FDR005

#Repeat for old only

#Subset for old
NSQ04_RKMP_R1_Old <- subset_samples(NSQ04_RKMP_R1_ni, SynCom.ID == "RKMP.Old")
sort(sample_sums(NSQ04_RKMP_R1_Old))
#extract count data 
NSQ04_RKMP_genotype_R1_Old <- otu_table(NSQ04_RKMP_R1_Old)
countData = as.data.frame(NSQ04_RKMP_genotype_R1_Old)
colnames(NSQ04_RKMP_genotype_R1_Old)
#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R1_Old)))
rownames(colData)
#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_RKMP_R1_Old_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)
#execute the differential count analysis with the function DESeq 
NSQ04_RKMP_R1_Old_test <- DESeq(NSQ04_RKMP_R1_Old_cds, fitType="local", betaPrior = FALSE) 
#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_R1_Old_Ba <- results(NSQ04_RKMP_R1_Old_test, contrast = c("Genotype", "Bulk", "Barke")) 
#extract  ASVs whose adjusted p.value in a given comparison is below 0.05 [if any]
rhizo_genotype_RKMP_R1_Old_Ba_FDR005 <- rhizo_genotype_RKMP_R1_Old_Ba[(rownames(rhizo_genotype_RKMP_R1_Old_Ba)[which(rhizo_genotype_RKMP_R1_Old_Ba$padj <0.05)]), ]
#inspect the generated file
rhizo_genotype_RKMP_R1_Old_Ba_FDR005

#Proportion charts

#set the colours: no inoculum, white bulk
SynCom_col_ni_wb <- c("#FFFFFF", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_R1_Old_Ba <- results(NSQ04_RKMP_R1_Old_test, contrast = c("Genotype", "Bulk", "Barke"))

#extract  first row
rhizo_genotype_r1 <- rhizo_genotype_RKMP_R1_Old_Ba[1, ]
#inspect the generated file
rhizo_genotype_r1
#Link back to full dataset
RKMP_r1 <- prune_taxa(rownames(rhizo_genotype_r1), NSQ04_RKMP_R1_Old)
RKMP_r1
#who is there?
tax_table(RKMP_r1)
#Pseudomonas
#% of enriched vs. total
r1_enriched <- as.data.frame(sample_sums(RKMP_r1)/sample_sums(NSQ04_RKMP_R1_Old))
colnames(r1_enriched) <- c("Pseudomonas_proportion")
#extract the mapping file
r1_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R1_Old)))
#merge the dataset
r1_plotting <- cbind(r1_map, r1_enriched)
#Order the factors
r1_plotting$Genotype <- ordered(r1_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
p <- ggplot(r1_plotting, aes(x=Genotype, y=Pseudomonas_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb)
p + geom_point(aes(y=Pseudomonas_proportion, group=Genotype), position = position_dodge(width=0.75))

#extract row 2
rhizo_genotype_r2 <- rhizo_genotype_RKMP_R1_Old_Ba[2, ]
#inspect the generated file
rhizo_genotype_r2
#Link back to full dataset
RKMP_r2 <- prune_taxa(rownames(rhizo_genotype_r2), NSQ04_RKMP_R1_Old)
RKMP_r2
#who is there?
tax_table(RKMP_r2)
#Pedobacter
#% of enriched vs. total
r2_enriched <- as.data.frame(sample_sums(RKMP_r2)/sample_sums(NSQ04_RKMP_R1_Old))
colnames(r2_enriched) <- c("Pedobacter_proportion")
#extract the mapping file
r2_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R1_Old)))
#merge the dataset
r2_plotting <- cbind(r2_map, r2_enriched)
#Order the factors
r2_plotting$Genotype <- ordered(r2_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
p <- ggplot(r2_plotting, aes(x=Genotype, y=Pedobacter_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb)
p + geom_point(aes(y=Pedobacter_proportion, group=Genotype), position = position_dodge(width=0.75))

#extract row 3
rhizo_genotype_r3 <- rhizo_genotype_RKMP_R1_Old_Ba[3, ]
#inspect the generated file
rhizo_genotype_r3
#Link back to full dataset
RKMP_r3 <- prune_taxa(rownames(rhizo_genotype_r3), NSQ04_RKMP_R1_Old)
RKMP_r3
#who is there?
tax_table(RKMP_r3)
#Chryseobacterium
#% of enriched vs. total
r3_enriched <- as.data.frame(sample_sums(RKMP_r3)/sample_sums(NSQ04_RKMP_R1_Old))
colnames(r3_enriched) <- c("Chryseobacterium_proportion")
#extract the mapping file
r3_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R1_Old)))
#merge the dataset
r3_plotting <- cbind(r3_map, r3_enriched)
#Order the factors
r3_plotting$Genotype <- ordered(r3_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
p <- ggplot(r3_plotting, aes(x=Genotype, y=Chryseobacterium_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb)
p + geom_point(aes(y=Chryseobacterium_proportion, group=Genotype), position = position_dodge(width=0.75))

#extract row 4
rhizo_genotype_r4 <- rhizo_genotype_RKMP_R1_Old_Ba[4, ]
#inspect the generated file
rhizo_genotype_r4
#Link back to full dataset
RKMP_r4 <- prune_taxa(rownames(rhizo_genotype_r4), NSQ04_RKMP_R1_Old)
RKMP_r4
#who is there?
tax_table(RKMP_r4)
#Stenotrophomonas
#% of enriched vs. total
r4_enriched <- as.data.frame(sample_sums(RKMP_r4)/sample_sums(NSQ04_RKMP_R1_Old))
colnames(r4_enriched) <- c("Stenotrophomonas_proportion")
#extract the mapping file
r4_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R1_Old)))
#merge the dataset
r4_plotting <- cbind(r4_map, r4_enriched)
#Order the factors
r4_plotting$Genotype <- ordered(r4_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
p <- ggplot(r4_plotting, aes(x=Genotype, y=Stenotrophomonas_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb)
p + geom_point(aes(y=Stenotrophomonas_proportion, group=Genotype), position = position_dodge(width=0.75))

#extract row 5
rhizo_genotype_r5 <- rhizo_genotype_RKMP_R1_Old_Ba[5, ]
#inspect the generated file
rhizo_genotype_r5
#Link back to full dataset
RKMP_r5 <- prune_taxa(rownames(rhizo_genotype_r5), NSQ04_RKMP_R1_Old)
RKMP_r5
#who is there?
tax_table(RKMP_r5)
#Bacillus
#% of enriched vs. total
r5_enriched <- as.data.frame(sample_sums(RKMP_r5)/sample_sums(NSQ04_RKMP_R1_Old))
colnames(r5_enriched) <- c("Bacillus_proportion")
#extract the mapping file
r5_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R1_Old)))
#merge the dataset
r5_plotting <- cbind(r5_map, r5_enriched)
#Order the factors
r5_plotting$Genotype <- ordered(r5_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
p <- ggplot(r5_plotting, aes(x=Genotype, y=Bacillus_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb)
p + geom_point(aes(y=Bacillus_proportion, group=Genotype), position = position_dodge(width=0.75))

#Repeat for New only

#Subset for New
NSQ04_RKMP_R1_New <- subset_samples(NSQ04_RKMP_R1_ni, SynCom.ID == "RKMP.New")
sort(sample_sums(NSQ04_RKMP_R1_New))
#extract count data 
NSQ04_RKMP_genotype_R1_New <- otu_table(NSQ04_RKMP_R1_New)
countData = as.data.frame(NSQ04_RKMP_genotype_R1_New)
colnames(NSQ04_RKMP_genotype_R1_New)
#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R1_New)))
rownames(colData)
#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_RKMP_R1_New_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)
#execute the differential count analysis with the function DESeq 
NSQ04_RKMP_R1_New_test <- DESeq(NSQ04_RKMP_R1_New_cds, fitType="local", betaPrior = FALSE) 
#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_R1_New_Ba <- results(NSQ04_RKMP_R1_New_test, contrast = c("Genotype", "Bulk", "Barke")) 
#extract  ASVs whose adjusted p.value in a given comparison is below 0.05 [if any]
rhizo_genotype_RKMP_R1_New_Ba_FDR005 <- rhizo_genotype_RKMP_R1_New_Ba[(rownames(rhizo_genotype_RKMP_R1_New_Ba)[which(rhizo_genotype_RKMP_R1_New_Ba$padj <0.05)]), ]
#inspect the generated file
rhizo_genotype_RKMP_R1_New_Ba_FDR005

#Proportion charts

#set the colours: no inoculum, white bulk
SynCom_col_ni_wb <- c("#FFFFFF", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_R1_New_Ba <- results(NSQ04_RKMP_R1_New_test, contrast = c("Genotype", "Bulk", "Barke"))

#extract  first row
rhizo_genotype_r1 <- rhizo_genotype_RKMP_R1_New_Ba[1, ]
#inspect the generated file
rhizo_genotype_r1
#Link back to full dataset
RKMP_r1 <- prune_taxa(rownames(rhizo_genotype_r1), NSQ04_RKMP_R1_New)
RKMP_r1
#who is there?
tax_table(RKMP_r1)
#Pseudomonas
#% of enriched vs. total
r1_enriched <- as.data.frame(sample_sums(RKMP_r1)/sample_sums(NSQ04_RKMP_R1_New))
colnames(r1_enriched) <- c("Pseudomonas_proportion")
#extract the mapping file
r1_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R1_New)))
#merge the dataset
r1_plotting <- cbind(r1_map, r1_enriched)
#Order the factors
r1_plotting$Genotype <- ordered(r1_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
p <- ggplot(r1_plotting, aes(x=Genotype, y=Pseudomonas_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb)
p + geom_point(aes(y=Pseudomonas_proportion, group=Genotype), position = position_dodge(width=0.75))

#extract row 2
rhizo_genotype_r2 <- rhizo_genotype_RKMP_R1_New_Ba[2, ]
#inspect the generated file
rhizo_genotype_r2
#Link back to full dataset
RKMP_r2 <- prune_taxa(rownames(rhizo_genotype_r2), NSQ04_RKMP_R1_New)
RKMP_r2
#who is there?
tax_table(RKMP_r2)
#Pedobacter
#% of enriched vs. total
r2_enriched <- as.data.frame(sample_sums(RKMP_r2)/sample_sums(NSQ04_RKMP_R1_New))
colnames(r2_enriched) <- c("Pedobacter_proportion")
#extract the mapping file
r2_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R1_New)))
#merge the dataset
r2_plotting <- cbind(r2_map, r2_enriched)
#Order the factors
r2_plotting$Genotype <- ordered(r2_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
p <- ggplot(r2_plotting, aes(x=Genotype, y=Pedobacter_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb)
p + geom_point(aes(y=Pedobacter_proportion, group=Genotype), position = position_dodge(width=0.75))

#extract row 3
rhizo_genotype_r3 <- rhizo_genotype_RKMP_R1_New_Ba[3, ]
#inspect the generated file
rhizo_genotype_r3
#Link back to full dataset
RKMP_r3 <- prune_taxa(rownames(rhizo_genotype_r3), NSQ04_RKMP_R1_New)
RKMP_r3
#who is there?
tax_table(RKMP_r3)
#Chryseobacterium
#% of enriched vs. total
r3_enriched <- as.data.frame(sample_sums(RKMP_r3)/sample_sums(NSQ04_RKMP_R1_New))
colnames(r3_enriched) <- c("Chryseobacterium_proportion")
#extract the mapping file
r3_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R1_New)))
#merge the dataset
r3_plotting <- cbind(r3_map, r3_enriched)
#Order the factors
r3_plotting$Genotype <- ordered(r3_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
p <- ggplot(r3_plotting, aes(x=Genotype, y=Chryseobacterium_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb)
p + geom_point(aes(y=Chryseobacterium_proportion, group=Genotype), position = position_dodge(width=0.75))

#extract row 4
rhizo_genotype_r4 <- rhizo_genotype_RKMP_R1_New_Ba[4, ]
#inspect the generated file
rhizo_genotype_r4
#Link back to full dataset
RKMP_r4 <- prune_taxa(rownames(rhizo_genotype_r4), NSQ04_RKMP_R1_New)
RKMP_r4
#who is there?
tax_table(RKMP_r4)
#Stenotrophomonas
#% of enriched vs. total
r4_enriched <- as.data.frame(sample_sums(RKMP_r4)/sample_sums(NSQ04_RKMP_R1_New))
colnames(r4_enriched) <- c("Stenotrophomonas_proportion")
#extract the mapping file
r4_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R1_New)))
#merge the dataset
r4_plotting <- cbind(r4_map, r4_enriched)
#Order the factors
r4_plotting$Genotype <- ordered(r4_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
p <- ggplot(r4_plotting, aes(x=Genotype, y=Stenotrophomonas_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb)
p + geom_point(aes(y=Stenotrophomonas_proportion, group=Genotype), position = position_dodge(width=0.75))

#extract row 5
rhizo_genotype_r5 <- rhizo_genotype_RKMP_R1_New_Ba[5, ]
#inspect the generated file
rhizo_genotype_r5
#Link back to full dataset
RKMP_r5 <- prune_taxa(rownames(rhizo_genotype_r5), NSQ04_RKMP_R1_New)
RKMP_r5
#who is there?
tax_table(RKMP_r5)
#Bacillus
#% of enriched vs. total
r5_enriched <- as.data.frame(sample_sums(RKMP_r5)/sample_sums(NSQ04_RKMP_R1_New))
colnames(r5_enriched) <- c("Bacillus_proportion")
#extract the mapping file
r5_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R1_New)))
#merge the dataset
r5_plotting <- cbind(r5_map, r5_enriched)
#Order the factors
r5_plotting$Genotype <- ordered(r5_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
p <- ggplot(r5_plotting, aes(x=Genotype, y=Bacillus_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb)
p + geom_point(aes(y=Bacillus_proportion, group=Genotype), position = position_dodge(width=0.75))

##############################################################
#Run 2 Only
##############################################################

#re-order the factors
sample_data(NSQ04_RKMP_calculation)$Genotype <- factor(sample_data(NSQ04_RKMP_calculation)$Genotype, levels=c("Inoculum", "Bulk", "Barke", "124-17", "124-52", "Morex"))
sample_data(NSQ04_RKMP_calculation)$SynCom.ID <- factor(sample_data(NSQ04_RKMP_calculation)$SynCom.ID, levels=c("RKMP.Old", "RKMP.New"))
sample_data(NSQ04_RKMP_calculation)$Plant <- factor(sample_data(NSQ04_RKMP_calculation)$Plant, levels=c("Unplanted", "Planted"))
sample_data(NSQ04_RKMP_calculation)$Date.Inoc <- as.factor(sample_data(NSQ04_RKMP_calculation)$Date.Inoc)

#set the colours: with inoculum
SynCom_col <- c("#999999", "#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")
#set the colours: no inoculum
SynCom_col_ni <- c("#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")
#set the colours: plant
Plant_col <- c("#E1BE6A", "#40B0A6")
#set the colours: no inoculum; no bulk
SynCom_col_plant <- c("#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#Subset for run 2
NSQ04_RKMP_R2 <- subset_samples(NSQ04_RKMP_calculation, Date.Inoc == "241010")
sort(sample_sums(NSQ04_RKMP_R2))

#Stacked bar plot
p = plot_bar(NSQ04_RKMP_R2, "Genus", facet_grid=SynCom.ID ~ Genotype, fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#Without inoculum
NSQ04_RKMP_R2_ni <- subset_samples(NSQ04_RKMP_R2, Genotype != "Inoculum")
sort(sample_sums(NSQ04_RKMP_R2_ni))

#Stacked bar plot
p = plot_bar(NSQ04_RKMP_R2_ni, "Genus", facet_grid=SynCom.ID ~ Genotype, fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#alpha diversity calculation
#https://www.nature.com/articles/s41598-024-77864-y
p = plot_richness(NSQ04_RKMP_R2, x="SynCom.ID", color="Genotype", measures=c("Observed", "Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#remove Observed, this is just a confirmation
p = plot_richness(NSQ04_RKMP_R2, x="SynCom.ID", color="Genotype", measures=c("Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#without inoculum
p = plot_richness(NSQ04_RKMP_R2_ni, x="SynCom.ID", color="Genotype", measures=c("Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col_ni)
p

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
NSQ04_RKMP_R2.ord <- ordinate(NSQ04_RKMP_R2, "NMDS", "bray")
#plot
p = plot_ordination(NSQ04_RKMP_R2, NSQ04_RKMP_R2.ord, type="samples", color="Genotype", shape="SynCom.ID") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#unconstrained ordination - no inoculum
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
NSQ04_RKMP_R2_ni.ord <- ordinate(NSQ04_RKMP_R2_ni, "NMDS", "bray")
#plot
p = plot_ordination(NSQ04_RKMP_R2_ni, NSQ04_RKMP_R2_ni.ord, type="samples", color="Genotype", shape="SynCom.ID") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col_ni)
p

#Planted only
NSQ04_RKMP_R2_plant <- subset_samples(NSQ04_RKMP_R2_ni, Plant == "Planted")
sort(sample_sums(NSQ04_RKMP_R2_plant))

#unconstrained ordination - planted only
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
NSQ04_RKMP_R2_plant.ord <- ordinate(NSQ04_RKMP_R2_plant, "NMDS", "bray")
#plot
p = plot_ordination(NSQ04_RKMP_R2_plant, NSQ04_RKMP_R2_plant.ord, type="samples", color="Genotype", shape="SynCom.ID") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col_plant)
p

#Constrained for SynCom effect
NSQ04_RKMP_R2_S.cap <- ordinate(NSQ04_RKMP_R2, "CAP", "bray", ~ SynCom.ID)
#Plot
a = plot_ordination(NSQ04_RKMP_R2, NSQ04_RKMP_R2_S.cap, type="samples", color="Genotype", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Constrained for SynCom effect - no inoculum
NSQ04_RKMP_R2_ni_S.cap <- ordinate(NSQ04_RKMP_R2_ni, "CAP", "bray", ~ SynCom.ID)
#Plot
a = plot_ordination(NSQ04_RKMP_R2_ni, NSQ04_RKMP_R2_ni_S.cap, type="samples", color="Genotype", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_ni)
a

#Constrained for SynCom effect - planted
NSQ04_RKMP_R2_plant_S.cap <- ordinate(NSQ04_RKMP_R2_plant, "CAP", "bray", ~ SynCom.ID)
#Plot
a = plot_ordination(NSQ04_RKMP_R2_plant, NSQ04_RKMP_R2_plant_S.cap, type="samples", color="Genotype", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_plant)
a

#Constrained for Genotype effect
NSQ04_RKMP_R2_G.cap <- ordinate(NSQ04_RKMP_R2, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(NSQ04_RKMP_R2, NSQ04_RKMP_R2_G.cap, type="samples", color="Genotype", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Constrained for genotype effect - no inoculum
NSQ04_RKMP_R2_ni_P.cap <- ordinate(NSQ04_RKMP_R2_ni, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(NSQ04_RKMP_R2_ni, NSQ04_RKMP_R2_ni_P.cap, type="samples", color="Genotype", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_ni)
a

#Stats - Permanova
#Plant Effect
BC_bacteria  <- phyloseq::distance(NSQ04_RKMP_R2_ni, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Plant * SynCom.ID, data= as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R2_ni))), permutations = 5000)
Stat
#Plant Effect Plus Genotype effect
BC_bacteria  <- phyloseq::distance(NSQ04_RKMP_R2_ni, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Plant * SynCom.ID * Genotype, data= as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R2_ni))), permutations = 5000)
Stat

####PAPER####
#Constrained for genotype effect - planted only
NSQ04_RKMP_R2_plant_P.cap <- ordinate(NSQ04_RKMP_R2_plant, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(NSQ04_RKMP_R2_plant, NSQ04_RKMP_R2_plant_P.cap, type="samples", color="Genotype", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_plant)+
  scale_shape_discrete(name = "Assembly Method", labels = c("By OD", "By CFU")) + theme(legend.position="bottom") 
a

#Stats - Permanova
#Genotype Effect
BC_bacteria  <- phyloseq::distance(NSQ04_RKMP_R2_plant, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Genotype * SynCom.ID, data= as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R2_plant))), permutations = 5000)
Stat

#Deseq
#extract count data 
NSQ04_RKMP_genotype_R2_ni <- otu_table(NSQ04_RKMP_R2_ni)
countData = as.data.frame(NSQ04_RKMP_genotype_R2_ni)
colnames(NSQ04_RKMP_genotype_R2_ni)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R2_ni)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_RKMP_R2_ni_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)

#execute the differential count analysis with the function DESeq 
NSQ04_RKMP_R2_ni_test <- DESeq(NSQ04_RKMP_R2_ni_cds, fitType="local", betaPrior = FALSE) 

#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_R2_ni_Ba <- results(NSQ04_RKMP_R2_ni_test, contrast = c("Genotype", "Bulk", "Barke")) 

#extract  ASVs whose adjusted p.value in a given comparison is below 0.05 [if any]
rhizo_genotype_RKMP_R2_ni_Ba_FDR005 <- rhizo_genotype_RKMP_R2_ni_Ba[(rownames(rhizo_genotype_RKMP_R2_ni_Ba)[which(rhizo_genotype_RKMP_R2_ni_Ba$padj <0.05)]), ]

#inspect the generated file
rhizo_genotype_RKMP_R2_ni_Ba_FDR005

#Repeat for 17
rhizo_genotype_RKMP_R2_ni_17 <- results(NSQ04_RKMP_R2_ni_test, contrast = c("Genotype", "Bulk", "124-17")) 
rhizo_genotype_RKMP_R2_ni_17_FDR005 <- rhizo_genotype_RKMP_R2_ni_17[(rownames(rhizo_genotype_RKMP_R2_ni_17)[which(rhizo_genotype_RKMP_R2_ni_17$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_ni_17_FDR005

#Repeat for 52
rhizo_genotype_RKMP_R2_ni_52 <- results(NSQ04_RKMP_R2_ni_test, contrast = c("Genotype", "Bulk", "124-52")) 
rhizo_genotype_RKMP_R2_ni_52_FDR005 <- rhizo_genotype_RKMP_R2_ni_52[(rownames(rhizo_genotype_RKMP_R2_ni_52)[which(rhizo_genotype_RKMP_R2_ni_52$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_ni_52_FDR005

#Repeat for Morex
rhizo_genotype_RKMP_R2_ni_Mx <- results(NSQ04_RKMP_R2_ni_test, contrast = c("Genotype", "Bulk", "Morex")) 
rhizo_genotype_RKMP_R2_ni_Mx_FDR005 <- rhizo_genotype_RKMP_R2_ni_Mx[(rownames(rhizo_genotype_RKMP_R2_ni_Mx)[which(rhizo_genotype_RKMP_R2_ni_Mx$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_ni_Mx_FDR005

#Switch to Barke as reference point

#Repeat for 17
rhizo_genotype_RKMP_R2_ni_17 <- results(NSQ04_RKMP_R2_ni_test, contrast = c("Genotype", "Barke", "124-17")) 
rhizo_genotype_RKMP_R2_ni_17_FDR005 <- rhizo_genotype_RKMP_R2_ni_17[(rownames(rhizo_genotype_RKMP_R2_ni_17)[which(rhizo_genotype_RKMP_R2_ni_17$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_ni_17_FDR005

#Repeat for 52
rhizo_genotype_RKMP_R2_ni_52 <- results(NSQ04_RKMP_R2_ni_test, contrast = c("Genotype", "Barke", "124-52")) 
rhizo_genotype_RKMP_R2_ni_52_FDR005 <- rhizo_genotype_RKMP_R2_ni_52[(rownames(rhizo_genotype_RKMP_R2_ni_52)[which(rhizo_genotype_RKMP_R2_ni_52$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_ni_52_FDR005

#Repeat for Morex
rhizo_genotype_RKMP_R2_ni_Mx <- results(NSQ04_RKMP_R2_ni_test, contrast = c("Genotype", "Barke", "Morex")) 
rhizo_genotype_RKMP_R2_ni_Mx_FDR005 <- rhizo_genotype_RKMP_R2_ni_Mx[(rownames(rhizo_genotype_RKMP_R2_ni_Mx)[which(rhizo_genotype_RKMP_R2_ni_Mx$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_ni_Mx_FDR005

#Repeat for old only

#Subset for old
NSQ04_RKMP_R2_Old <- subset_samples(NSQ04_RKMP_R2_ni, SynCom.ID == "RKMP.Old")
sort(sample_sums(NSQ04_RKMP_R2_Old))
#extract count data 
NSQ04_RKMP_genotype_R2_Old <- otu_table(NSQ04_RKMP_R2_Old)
countData = as.data.frame(NSQ04_RKMP_genotype_R2_Old)
colnames(NSQ04_RKMP_genotype_R2_Old)
#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R2_Old)))
rownames(colData)
#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_RKMP_R2_Old_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)
#execute the differential count analysis with the function DESeq 
NSQ04_RKMP_R2_Old_test <- DESeq(NSQ04_RKMP_R2_Old_cds, fitType="local", betaPrior = FALSE) 

#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_R2_Old_Ba <- results(NSQ04_RKMP_R2_Old_test, contrast = c("Genotype", "Bulk", "Barke")) 
#extract  ASVs whose adjusted p.value in a given comparison is below 0.05 [if any]
rhizo_genotype_RKMP_R2_Old_Ba_FDR005 <- rhizo_genotype_RKMP_R2_Old_Ba[(rownames(rhizo_genotype_RKMP_R2_Old_Ba)[which(rhizo_genotype_RKMP_R2_Old_Ba$padj <0.05)]), ]
#inspect the generated file
rhizo_genotype_RKMP_R2_Old_Ba_FDR005

#17
rhizo_genotype_RKMP_R2_Old_17 <- results(NSQ04_RKMP_R2_Old_test, contrast = c("Genotype", "Bulk", "124-17")) 
rhizo_genotype_RKMP_R2_Old_17_FDR005 <- rhizo_genotype_RKMP_R2_Old_17[(rownames(rhizo_genotype_RKMP_R2_Old_17)[which(rhizo_genotype_RKMP_R2_Old_17$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_Old_17_FDR005

#52
rhizo_genotype_RKMP_R2_Old_52 <- results(NSQ04_RKMP_R2_Old_test, contrast = c("Genotype", "Bulk", "124-52")) 
rhizo_genotype_RKMP_R2_Old_52_FDR005 <- rhizo_genotype_RKMP_R2_Old_52[(rownames(rhizo_genotype_RKMP_R2_Old_52)[which(rhizo_genotype_RKMP_R2_Old_52$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_Old_52_FDR005

#Mx
rhizo_genotype_RKMP_R2_Old_Mx <- results(NSQ04_RKMP_R2_Old_test, contrast = c("Genotype", "Bulk", "Morex")) 
rhizo_genotype_RKMP_R2_Old_Mx_FDR005 <- rhizo_genotype_RKMP_R2_Old_Mx[(rownames(rhizo_genotype_RKMP_R2_Old_Mx)[which(rhizo_genotype_RKMP_R2_Old_Mx$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_Old_Mx_FDR005

#Switch to Barke as reference point

#Repeat for 17
rhizo_genotype_RKMP_R2_Old_17 <- results(NSQ04_RKMP_R2_Old_test, contrast = c("Genotype", "Barke", "124-17")) 
rhizo_genotype_RKMP_R2_Old_17_FDR005 <- rhizo_genotype_RKMP_R2_Old_17[(rownames(rhizo_genotype_RKMP_R2_Old_17)[which(rhizo_genotype_RKMP_R2_Old_17$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_Old_17_FDR005

#Repeat for 52
rhizo_genotype_RKMP_R2_Old_52 <- results(NSQ04_RKMP_R2_Old_test, contrast = c("Genotype", "Barke", "124-52")) 
rhizo_genotype_RKMP_R2_Old_52_FDR005 <- rhizo_genotype_RKMP_R2_Old_52[(rownames(rhizo_genotype_RKMP_R2_Old_52)[which(rhizo_genotype_RKMP_R2_Old_52$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_Old_52_FDR005

#Repeat for Morex
rhizo_genotype_RKMP_R2_Old_Mx <- results(NSQ04_RKMP_R2_Old_test, contrast = c("Genotype", "Barke", "Morex")) 
rhizo_genotype_RKMP_R2_Old_Mx_FDR005 <- rhizo_genotype_RKMP_R2_Old_Mx[(rownames(rhizo_genotype_RKMP_R2_Old_Mx)[which(rhizo_genotype_RKMP_R2_Old_Mx$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_Old_Mx_FDR005

####PAPER####
#Proportion charts

#set the colours: no inoculum, white bulk
SynCom_col_ni_wb <- c("#FFFFFF", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_R2_Old_Ba <- results(NSQ04_RKMP_R2_Old_test, contrast = c("Genotype", "Bulk", "Barke"))

#extract  first row
rhizo_genotype_r1 <- rhizo_genotype_RKMP_R2_Old_Ba[1, ]
#inspect the generated file
rhizo_genotype_r1
#Link back to full dataset
RKMP_r1 <- prune_taxa(rownames(rhizo_genotype_r1), NSQ04_RKMP_R2_Old)
RKMP_r1
#who is there?
tax_table(RKMP_r1)
#Pseudomonas
#% of enriched vs. total
r1_enriched <- as.data.frame(sample_sums(RKMP_r1)/sample_sums(NSQ04_RKMP_R2_Old))
colnames(r1_enriched) <- c("Pseudomonas_proportion")
#extract the mapping file
r1_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R2_Old)))
#merge the dataset
r1_plotting <- cbind(r1_map, r1_enriched)
#Order the factors
r1_plotting$Genotype <- ordered(r1_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
r1 <- ggplot(r1_plotting, aes(x=Genotype, y=Pseudomonas_proportion, fill = Genotype)) 
r1 <- r1 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb) + 
  geom_point(aes(y=Pseudomonas_proportion, group=Genotype), position = position_dodge(width=0.75))+
  ylim(0,0.8) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r1 

#extract row 2
rhizo_genotype_r2 <- rhizo_genotype_RKMP_R2_Old_Ba[2, ]
#inspect the generated file
rhizo_genotype_r2
#Link back to full dataset
RKMP_r2 <- prune_taxa(rownames(rhizo_genotype_r2), NSQ04_RKMP_R2_Old)
RKMP_r2
#who is there?
tax_table(RKMP_r2)
#Pedobacter
#% of enriched vs. total
r2_enriched <- as.data.frame(sample_sums(RKMP_r2)/sample_sums(NSQ04_RKMP_R2_Old))
colnames(r2_enriched) <- c("Pedobacter_proportion")
#extract the mapping file
r2_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R2_Old)))
#merge the dataset
r2_plotting <- cbind(r2_map, r2_enriched)
#Order the factors
r2_plotting$Genotype <- ordered(r2_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
r2 <- ggplot(r2_plotting, aes(x=Genotype, y=Pedobacter_proportion, fill = Genotype)) 
r2 <- r2 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb) + 
  geom_point(aes(y=Pedobacter_proportion, group=Genotype), position = position_dodge(width=0.75))+
  ylim(0,0.8) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r2 

#extract row 3
rhizo_genotype_r3 <- rhizo_genotype_RKMP_R2_Old_Ba[3, ]
#inspect the generated file
rhizo_genotype_r3
#Link back to full dataset
RKMP_r3 <- prune_taxa(rownames(rhizo_genotype_r3), NSQ04_RKMP_R2_Old)
RKMP_r3
#who is there?
tax_table(RKMP_r3)
#Chryseobacterium
#% of enriched vs. total
r3_enriched <- as.data.frame(sample_sums(RKMP_r3)/sample_sums(NSQ04_RKMP_R2_Old))
colnames(r3_enriched) <- c("Chryseobacterium_proportion")
#extract the mapping file
r3_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R2_Old)))
#merge the dataset
r3_plotting <- cbind(r3_map, r3_enriched)
#Order the factors
r3_plotting$Genotype <- ordered(r3_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
r3 <- ggplot(r3_plotting, aes(x=Genotype, y=Chryseobacterium_proportion, fill = Genotype)) 
r3 <- r3 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb) + 
  geom_point(aes(y=Chryseobacterium_proportion, group=Genotype), position = position_dodge(width=0.75)) +
  ylim(0,0.8) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r3 

#extract row 4
rhizo_genotype_r4 <- rhizo_genotype_RKMP_R2_Old_Ba[4, ]
#inspect the generated file
rhizo_genotype_r4
#Link back to full dataset
RKMP_r4 <- prune_taxa(rownames(rhizo_genotype_r4), NSQ04_RKMP_R2_Old)
RKMP_r4
#who is there?
tax_table(RKMP_r4)
#Stenotrophomonas
#% of enriched vs. total
r4_enriched <- as.data.frame(sample_sums(RKMP_r4)/sample_sums(NSQ04_RKMP_R2_Old))
colnames(r4_enriched) <- c("Stenotrophomonas_proportion")
#extract the mapping file
r4_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R2_Old)))
#merge the dataset
r4_plotting <- cbind(r4_map, r4_enriched)
#Order the factors
r4_plotting$Genotype <- ordered(r4_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
r4 <- ggplot(r4_plotting, aes(x=Genotype, y=Stenotrophomonas_proportion, fill = Genotype)) 
r4 <- r4 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb)+ 
  geom_point(aes(y=Stenotrophomonas_proportion, group=Genotype), position = position_dodge(width=0.75))+
  ylim(0,0.8) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r4 

#extract row 5
rhizo_genotype_r5 <- rhizo_genotype_RKMP_R2_Old_Ba[5, ]
#inspect the generated file
rhizo_genotype_r5
#Link back to full dataset
RKMP_r5 <- prune_taxa(rownames(rhizo_genotype_r5), NSQ04_RKMP_R2_Old)
RKMP_r5
#who is there?
tax_table(RKMP_r5)
#Bacillus
#% of enriched vs. total
r5_enriched <- as.data.frame(sample_sums(RKMP_r5)/sample_sums(NSQ04_RKMP_R2_Old))
colnames(r5_enriched) <- c("Bacillus_proportion")
#extract the mapping file
r5_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R2_Old)))
#merge the dataset
r5_plotting <- cbind(r5_map, r5_enriched)
#Order the factors
r5_plotting$Genotype <- ordered(r5_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
r5 <- ggplot(r5_plotting, aes(x=Genotype, y=Bacillus_proportion, fill = Genotype)) 
r5 <- r5 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb) + 
  geom_point(aes(y=Bacillus_proportion, group=Genotype), position = position_dodge(width=0.75))+
  ylim(0,0.8) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r5

#faceted graph of all
grid.arrange(r1, r2, r3, r4, r5, ncol=5)

#Repeat for New only

#Subset for New
NSQ04_RKMP_R2_New <- subset_samples(NSQ04_RKMP_R2_ni, SynCom.ID == "RKMP.New")
sort(sample_sums(NSQ04_RKMP_R2_New))
#extract count data 
NSQ04_RKMP_genotype_R2_New <- otu_table(NSQ04_RKMP_R2_New)
countData = as.data.frame(NSQ04_RKMP_genotype_R2_New)
colnames(NSQ04_RKMP_genotype_R2_New)
#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R2_New)))
rownames(colData)
#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_RKMP_R2_New_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)
#execute the differential count analysis with the function DESeq 
NSQ04_RKMP_R2_New_test <- DESeq(NSQ04_RKMP_R2_New_cds, fitType="local", betaPrior = FALSE) 

#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_R2_New_Ba <- results(NSQ04_RKMP_R2_New_test, contrast = c("Genotype", "Bulk", "Barke")) 
#extract  ASVs whose adjusted p.value in a given comparison is below 0.05 [if any]
rhizo_genotype_RKMP_R2_New_Ba_FDR005 <- rhizo_genotype_RKMP_R2_New_Ba[(rownames(rhizo_genotype_RKMP_R2_New_Ba)[which(rhizo_genotype_RKMP_R2_New_Ba$padj <0.05)]), ]
#inspect the generated file
rhizo_genotype_RKMP_R2_New_Ba_FDR005

#17
rhizo_genotype_RKMP_R2_New_17 <- results(NSQ04_RKMP_R2_New_test, contrast = c("Genotype", "Bulk", "124-17")) 
rhizo_genotype_RKMP_R2_New_17_FDR005 <- rhizo_genotype_RKMP_R2_New_17[(rownames(rhizo_genotype_RKMP_R2_New_17)[which(rhizo_genotype_RKMP_R2_New_17$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_New_17_FDR005

#52
rhizo_genotype_RKMP_R2_New_52 <- results(NSQ04_RKMP_R2_New_test, contrast = c("Genotype", "Bulk", "124-52")) 
rhizo_genotype_RKMP_R2_New_52_FDR005 <- rhizo_genotype_RKMP_R2_New_52[(rownames(rhizo_genotype_RKMP_R2_New_52)[which(rhizo_genotype_RKMP_R2_New_52$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_New_52_FDR005

#Mx
rhizo_genotype_RKMP_R2_New_Mx <- results(NSQ04_RKMP_R2_New_test, contrast = c("Genotype", "Bulk", "Morex")) 
rhizo_genotype_RKMP_R2_New_Mx_FDR005 <- rhizo_genotype_RKMP_R2_New_Mx[(rownames(rhizo_genotype_RKMP_R2_New_Mx)[which(rhizo_genotype_RKMP_R2_New_Mx$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_New_Mx_FDR005

#Switch to Barke as reference point

#Repeat for 17
rhizo_genotype_RKMP_R2_New_17 <- results(NSQ04_RKMP_R2_New_test, contrast = c("Genotype", "Barke", "124-17")) 
rhizo_genotype_RKMP_R2_New_17_FDR005 <- rhizo_genotype_RKMP_R2_New_17[(rownames(rhizo_genotype_RKMP_R2_New_17)[which(rhizo_genotype_RKMP_R2_New_17$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_New_17_FDR005

#Repeat for 52
rhizo_genotype_RKMP_R2_New_52 <- results(NSQ04_RKMP_R2_New_test, contrast = c("Genotype", "Barke", "124-52")) 
rhizo_genotype_RKMP_R2_New_52_FDR005 <- rhizo_genotype_RKMP_R2_New_52[(rownames(rhizo_genotype_RKMP_R2_New_52)[which(rhizo_genotype_RKMP_R2_New_52$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_New_52_FDR005

#Repeat for Morex
rhizo_genotype_RKMP_R2_New_Mx <- results(NSQ04_RKMP_R2_New_test, contrast = c("Genotype", "Barke", "Morex")) 
rhizo_genotype_RKMP_R2_New_Mx_FDR005 <- rhizo_genotype_RKMP_R2_New_Mx[(rownames(rhizo_genotype_RKMP_R2_New_Mx)[which(rhizo_genotype_RKMP_R2_New_Mx$padj <0.05)]), ]
rhizo_genotype_RKMP_R2_New_Mx_FDR005

####PAPER####
#Proportion charts

#set the colours: no inoculum, white bulk
SynCom_col_ni_wb <- c("#FFFFFF", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_R2_New_Ba <- results(NSQ04_RKMP_R2_New_test, contrast = c("Genotype", "Bulk", "Barke"))

#extract  first row
rhizo_genotype_r1 <- rhizo_genotype_RKMP_R2_New_Ba[1, ]
#inspect the generated file
rhizo_genotype_r1
#Link back to full dataset
RKMP_r1 <- prune_taxa(rownames(rhizo_genotype_r1), NSQ04_RKMP_R2_New)
RKMP_r1
#who is there?
tax_table(RKMP_r1)
#Pseudomonas
#% of enriched vs. total
r1_enriched <- as.data.frame(sample_sums(RKMP_r1)/sample_sums(NSQ04_RKMP_R2_New))
colnames(r1_enriched) <- c("Pseudomonas_proportion")
#extract the mapping file
r1_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R2_New)))
#merge the dataset
r1_plotting <- cbind(r1_map, r1_enriched)
#Order the factors
r1_plotting$Genotype <- ordered(r1_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
r1 <- ggplot(r1_plotting, aes(x=Genotype, y=Pseudomonas_proportion, fill = Genotype)) 
r1 <- r1 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb) + 
  geom_point(aes(y=Pseudomonas_proportion, group=Genotype), position = position_dodge(width=0.75))+
  ylim(0,0.8) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r1 

#extract row 2
rhizo_genotype_r2 <- rhizo_genotype_RKMP_R2_New_Ba[2, ]
#inspect the generated file
rhizo_genotype_r2
#Link back to full dataset
RKMP_r2 <- prune_taxa(rownames(rhizo_genotype_r2), NSQ04_RKMP_R2_New)
RKMP_r2
#who is there?
tax_table(RKMP_r2)
#Pedobacter
#% of enriched vs. total
r2_enriched <- as.data.frame(sample_sums(RKMP_r2)/sample_sums(NSQ04_RKMP_R2_New))
colnames(r2_enriched) <- c("Pedobacter_proportion")
#extract the mapping file
r2_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R2_New)))
#merge the dataset
r2_plotting <- cbind(r2_map, r2_enriched)
#Order the factors
r2_plotting$Genotype <- ordered(r2_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
r2 <- ggplot(r2_plotting, aes(x=Genotype, y=Pedobacter_proportion, fill = Genotype)) 
r2 <- r2 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb) + 
  geom_point(aes(y=Pedobacter_proportion, group=Genotype), position = position_dodge(width=0.75))+
  ylim(0,0.8) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r2 

#extract row 3
rhizo_genotype_r3 <- rhizo_genotype_RKMP_R2_New_Ba[3, ]
#inspect the generated file
rhizo_genotype_r3
#Link back to full dataset
RKMP_r3 <- prune_taxa(rownames(rhizo_genotype_r3), NSQ04_RKMP_R2_New)
RKMP_r3
#who is there?
tax_table(RKMP_r3)
#Chryseobacterium
#% of enriched vs. total
r3_enriched <- as.data.frame(sample_sums(RKMP_r3)/sample_sums(NSQ04_RKMP_R2_New))
colnames(r3_enriched) <- c("Chryseobacterium_proportion")
#extract the mapping file
r3_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R2_New)))
#merge the dataset
r3_plotting <- cbind(r3_map, r3_enriched)
#Order the factors
r3_plotting$Genotype <- ordered(r3_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
r3 <- ggplot(r3_plotting, aes(x=Genotype, y=Chryseobacterium_proportion, fill = Genotype)) 
r3 <- r3 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb) + 
  geom_point(aes(y=Chryseobacterium_proportion, group=Genotype), position = position_dodge(width=0.75))+
  ylim(0,0.8) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r3 

#extract row 4
rhizo_genotype_r4 <- rhizo_genotype_RKMP_R2_New_Ba[4, ]
#inspect the generated file
rhizo_genotype_r4
#Link back to full dataset
RKMP_r4 <- prune_taxa(rownames(rhizo_genotype_r4), NSQ04_RKMP_R2_New)
RKMP_r4
#who is there?
tax_table(RKMP_r4)
#Stenotrophomonas
#% of enriched vs. total
r4_enriched <- as.data.frame(sample_sums(RKMP_r4)/sample_sums(NSQ04_RKMP_R2_New))
colnames(r4_enriched) <- c("Stenotrophomonas_proportion")
#extract the mapping file
r4_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R2_New)))
#merge the dataset
r4_plotting <- cbind(r4_map, r4_enriched)
#Order the factors
r4_plotting$Genotype <- ordered(r4_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
r4 <- ggplot(r4_plotting, aes(x=Genotype, y=Stenotrophomonas_proportion, fill = Genotype)) 
r4 <- r4 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb)+ 
  geom_point(aes(y=Stenotrophomonas_proportion, group=Genotype), position = position_dodge(width=0.75))+
  ylim(0,0.8) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r4 

#extract row 5
rhizo_genotype_r5 <- rhizo_genotype_RKMP_R2_New_Ba[5, ]
#inspect the generated file
rhizo_genotype_r5
#Link back to full dataset
RKMP_r5 <- prune_taxa(rownames(rhizo_genotype_r5), NSQ04_RKMP_R2_New)
RKMP_r5
#who is there?
tax_table(RKMP_r5)
#Bacillus
#% of enriched vs. total
r5_enriched <- as.data.frame(sample_sums(RKMP_r5)/sample_sums(NSQ04_RKMP_R2_New))
colnames(r5_enriched) <- c("Bacillus_proportion")
#extract the mapping file
r5_map <- as.data.frame(as.matrix(sample_data(NSQ04_RKMP_R2_New)))
#merge the dataset
r5_plotting <- cbind(r5_map, r5_enriched)
#Order the factors
r5_plotting$Genotype <- ordered(r5_plotting$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
r5 <- ggplot(r5_plotting, aes(x=Genotype, y=Bacillus_proportion, fill = Genotype)) 
r5 <- r5 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb) + 
  geom_point(aes(y=Bacillus_proportion, group=Genotype), position = position_dodge(width=0.75))+
  ylim(0,0.8) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r5 

#faceted graph of all
grid.arrange(r1, r2, r3, r4, r5, ncol=5)

##############################################################
#Barke/Bulk Only
##############################################################

#re-order the factors
sample_data(NSQ04_RKMP_calculation)$Genotype <- factor(sample_data(NSQ04_RKMP_calculation)$Genotype, levels=c("Inoculum", "Bulk", "Barke", "124-17", "124-52", "Morex"))
sample_data(NSQ04_RKMP_calculation)$SynCom.ID <- factor(sample_data(NSQ04_RKMP_calculation)$SynCom.ID, levels=c("RKMP.Old", "RKMP.New"))
sample_data(NSQ04_RKMP_calculation)$Plant <- factor(sample_data(NSQ04_RKMP_calculation)$Plant, levels=c("Unplanted", "Planted"))
sample_data(NSQ04_RKMP_calculation)$Date.Inoc <- as.factor(sample_data(NSQ04_RKMP_calculation)$Date.Inoc)

#set the colours: with inoculum
SynCom_col <- c("#999999", "#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")
#set the colours: no inoculum
SynCom_col_ni <- c("#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")
#set the colours: plant
Plant_col <- c("#E1BE6A", "#40B0A6")

#Subset for Barke/Bulk
NSQ04_RKMP_BB <- subset_samples(NSQ04_RKMP_calculation, Genotype != "124-17")
NSQ04_RKMP_BB <- subset_samples(NSQ04_RKMP_BB, Genotype != "124-52")
NSQ04_RKMP_BB <- subset_samples(NSQ04_RKMP_BB, Genotype != "Morex")
sort(sample_sums(NSQ04_RKMP_BB))

#Subset for SynCom ID
NSQ04_RKMP_BB_Old <- subset_samples(NSQ04_RKMP_BB, SynCom.ID == "RKMP.Old")
NSQ04_RKMP_BB_New <- subset_samples(NSQ04_RKMP_BB, SynCom.ID == "RKMP.New")

#Stacked bar plot - old
o = plot_bar(NSQ04_RKMP_BB_Old, "Genus", facet_grid=Date.Inoc ~ Genotype, fill='Genus')
o = o + scale_fill_viridis(discrete=TRUE)
o = o + ggtitle('Old')
o

#Stacked bar plot - new
n = plot_bar(NSQ04_RKMP_BB_New, "Genus", facet_grid=Date.Inoc ~ Genotype, fill='Genus')
n = n + scale_fill_viridis(discrete=TRUE)
n = n + ggtitle('New')
n

#faceted graph of both
grid.arrange(o, n, ncol=2)

#Without inoculum
NSQ04_RKMP_BB_ni <- subset_samples(NSQ04_RKMP_BB, Genotype != "Inoculum")
NSQ04_RKMP_BB_ni_Old <- subset_samples(NSQ04_RKMP_BB_ni, SynCom.ID == "RKMP.Old")
NSQ04_RKMP_BB_ni_New <- subset_samples(NSQ04_RKMP_BB_ni, SynCom.ID == "RKMP.New")

#Stacked bar plot - old
p = plot_bar(NSQ04_RKMP_BB_ni_Old, "Genus", facet_grid=Date.Inoc ~ Genotype, fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#Stacked bar plot - new
p = plot_bar(NSQ04_RKMP_BB_ni_New, "Genus", facet_grid=Date.Inoc ~ Genotype, fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#alpha diversity calculation - Old
p = plot_richness(NSQ04_RKMP_BB_Old, x="Date.Inoc", color="Genotype", measures=c("Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#alpha diversity calculation - New
p = plot_richness(NSQ04_RKMP_BB_New, x="Date.Inoc", color="Genotype", measures=c("Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#without inoculum
#alpha diversity calculation - Old
p = plot_richness(NSQ04_RKMP_BB_ni_Old, x="Date.Inoc", color="Plant", measures=c("Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = Plant_col)
p

#alpha diversity calculation - New
p = plot_richness(NSQ04_RKMP_BB_ni_New, x="Date.Inoc", color="Plant", measures=c("Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = Plant_col)
p

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
NSQ04_RKMP_BB.ord <- ordinate(NSQ04_RKMP_BB, "NMDS", "bray")
#plot
p = plot_ordination(NSQ04_RKMP_BB, NSQ04_RKMP_BB.ord, type="samples", color="Genotype", shape="SynCom.ID") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p = p + facet_wrap(~Date.Inoc, 1)
p

#unconstrained ordination - no inoculum
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
NSQ04_RKMP_BB_ni.ord <- ordinate(NSQ04_RKMP_BB_ni, "NMDS", "bray")
#plot
p = plot_ordination(NSQ04_RKMP_BB_ni, NSQ04_RKMP_BB_ni.ord, type="samples", color="Plant", shape="SynCom.ID") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = Plant_col)
p = p + facet_wrap(~Date.Inoc, 1)
p

####PAPER####
#Constrained for SynCom effect
NSQ04_RKMP_BB_S.cap <- ordinate(NSQ04_RKMP_BB, "CAP", "bray", ~ SynCom.ID)
#Plot
a = plot_ordination(NSQ04_RKMP_BB, NSQ04_RKMP_BB_S.cap, type="samples", color="Genotype", shape="SynCom.ID") + 
  geom_point(size = 4, alpha = 0.75)+ 
  scale_shape_manual(values = c(16, 17))+
  scale_colour_manual(values = SynCom_col)+
  facet_wrap(~Date.Inoc, 1, labeller = labeller(Date.Inoc = c("241003" = "Replicate 1", "241010" = "Replicate 2")))+
  scale_shape_discrete(name = "Assembly Method", labels = c("By OD", "By CFU"))+ theme(legend.position="bottom")
a

#Constrained for SynCom effect - no inoculum
NSQ04_RKMP_BB_ni_S.cap <- ordinate(NSQ04_RKMP_BB_ni, "CAP", "bray", ~ SynCom.ID)
#Plot
a = plot_ordination(NSQ04_RKMP_BB_ni, NSQ04_RKMP_BB_ni_S.cap, type="samples", color="Plant", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = Plant_col)
a = a + facet_wrap(~Date.Inoc, 1)
a

#Constrained for 'Genotype' effect
NSQ04_RKMP_BB_G.cap <- ordinate(NSQ04_RKMP_BB, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(NSQ04_RKMP_BB, NSQ04_RKMP_BB_G.cap, type="samples", color="Genotype", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a = a + facet_wrap(~Date.Inoc, 1)
a

#Constrained for plant effect - no inoculum
NSQ04_RKMP_BB_ni_P.cap <- ordinate(NSQ04_RKMP_BB_ni, "CAP", "bray", ~ Plant)
#Plot
a = plot_ordination(NSQ04_RKMP_BB_ni, NSQ04_RKMP_BB_ni_P.cap, type="samples", color="Plant", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = Plant_col)
a = a + facet_wrap(~Date.Inoc, 1)
a

#Stats - Permanova
BC_bacteria  <- phyloseq::distance(NSQ04_RKMP_BB_ni, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Plant * SynCom.ID * Date.Inoc, data= as.data.frame(as.matrix(sample_data(NSQ04_RKMP_BB_ni))), permutations = 5000)
Stat

#Deseq
#extract count data 
NSQ04_RKMP_genotype_BB_ni <- otu_table(NSQ04_RKMP_BB_ni)
countData = as.data.frame(NSQ04_RKMP_genotype_BB_ni)
colnames(NSQ04_RKMP_genotype_BB_ni)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_RKMP_BB_ni)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_RKMP_BB_ni_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)

#execute the differential count analysis with the function DESeq 
NSQ04_RKMP_BB_ni_test <- DESeq(NSQ04_RKMP_BB_ni_cds, fitType="local", betaPrior = FALSE) 

#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_BB_ni_Ba <- results(NSQ04_RKMP_BB_ni_test, contrast = c("Genotype", "Bulk", "Barke")) 

#extract  ASVs whose adjusted p.value in a given comparison is below 0.05 [if any]
rhizo_genotype_RKMP_BB_ni_Ba_FDR005 <- rhizo_genotype_RKMP_BB_ni_Ba[(rownames(rhizo_genotype_RKMP_BB_ni_Ba)[which(rhizo_genotype_RKMP_BB_ni_Ba$padj <0.05)]), ]

#inspect the generated file
rhizo_genotype_RKMP_BB_ni_Ba_FDR005

#Repeat for old only

#extract count data 
NSQ04_RKMP_genotype_BB_Old <- otu_table(NSQ04_RKMP_BB_Old)
countData = as.data.frame(NSQ04_RKMP_genotype_BB_Old)
colnames(NSQ04_RKMP_genotype_BB_Old)
#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_RKMP_BB_Old)))
rownames(colData)
#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_RKMP_BB_Old_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)
#execute the differential count analysis with the function DESeq 
NSQ04_RKMP_BB_Old_test <- DESeq(NSQ04_RKMP_BB_Old_cds, fitType="local", betaPrior = FALSE) 
#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_BB_Old_Ba <- results(NSQ04_RKMP_BB_Old_test, contrast = c("Genotype", "Bulk", "Barke")) 
#extract  ASVs whose adjusted p.value in a given comparison is below 0.05 [if any]
rhizo_genotype_RKMP_BB_Old_Ba_FDR005 <- rhizo_genotype_RKMP_BB_Old_Ba[(rownames(rhizo_genotype_RKMP_BB_Old_Ba)[which(rhizo_genotype_RKMP_BB_Old_Ba$padj <0.05)]), ]
#inspect the generated file
rhizo_genotype_RKMP_BB_Old_Ba_FDR005

#Repeat for New only

#extract count data 
NSQ04_RKMP_genotype_BB_New <- otu_table(NSQ04_RKMP_BB_New)
countData = as.data.frame(NSQ04_RKMP_genotype_BB_New)
colnames(NSQ04_RKMP_genotype_BB_New)
#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_RKMP_BB_New)))
rownames(colData)
#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_RKMP_BB_New_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)
#execute the differential count analysis with the function DESeq 
NSQ04_RKMP_BB_New_test <- DESeq(NSQ04_RKMP_BB_New_cds, fitType="local", betaPrior = FALSE) 
#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_BB_New_Ba <- results(NSQ04_RKMP_BB_New_test, contrast = c("Genotype", "Bulk", "Barke")) 
#extract  ASVs whose adjusted p.value in a given comparison is below 0.05 [if any]
rhizo_genotype_RKMP_BB_New_Ba_FDR005 <- rhizo_genotype_RKMP_BB_New_Ba[(rownames(rhizo_genotype_RKMP_BB_New_Ba)[which(rhizo_genotype_RKMP_BB_New_Ba$padj <0.05)]), ]
#inspect the generated file
rhizo_genotype_RKMP_BB_New_Ba_FDR005


#########################################################################
# work with HK samples and set criteria to work with unevenly sequenced samples
# a) remove samples with less than 1,000 reads
# b) keep ASV with more than 15 reads in at least 10% of the samples
# c) rarefy at even sequencing depth
# d) agglomerate at genus level
#########################################################################

#Inspect
sort(sample_sums(NSQ04_RKMP_HK))

#remove poorly sequenced samples
NSQ04_RKMP_HK_1K <- prune_samples(sample_sums(NSQ04_RKMP_HK) > 1000, NSQ04_RKMP_HK)

#set an arbitrary threshold of 10 reads in 10% of samples
NSQ04_RKMP_HK_1K_threshold = filter_taxa(NSQ04_RKMP_HK_1K, function(x) sum(x > 10) > (0.1 *length(x)), TRUE)
NSQ04_RKMP_HK_1K_threshold

#inspect read distribution upon thresholding
sort(sample_sums(NSQ04_RKMP_HK_1K_threshold))

#agglomerate at genus level
NSQ04_RKMP_HK_genus <- tax_glom(NSQ04_RKMP_HK_1K_threshold, taxrank="Genus")

#inspect the files
NSQ04_RKMP_HK_genus
sort(sample_sums(NSQ04_RKMP_HK_genus))

#rarefy at 600
NSQ04_RKMP_HK_threshold_600 <- rarefy_even_depth(NSQ04_RKMP_HK_genus, 600)

#data visualisation
p = plot_bar(NSQ04_RKMP_HK_threshold_600, "Genus", facet_grid=SynCom.ID ~ Genotype, fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#save the file for the reproducibility of the code
#saveRDS(NSQ04_RKMP_HK_threshold_600, file = "NSQ04_RKMP_HK_calculation.rds")
NSQ04_RKMP_calculation_HK <-readRDS("NSQ04_RKMP_HK_calculation.rds")
NSQ04_RKMP_calculation_HK

#########
#HK Only#
#########

#re-order the factors
sample_data(NSQ04_RKMP_calculation_HK)$Genotype <- factor(sample_data(NSQ04_RKMP_calculation_HK)$Genotype, levels=c("Inoculum", "Bulk", "Barke", "124-17", "124-52", "Morex"))
sample_data(NSQ04_RKMP_calculation_HK)$SynCom.ID <- factor(sample_data(NSQ04_RKMP_calculation_HK)$SynCom.ID, levels=c("RKMP.Old", "RKMP.New"))
sample_data(NSQ04_RKMP_calculation_HK)$Plant <- factor(sample_data(NSQ04_RKMP_calculation_HK)$Plant, levels=c("Unplanted", "Planted"))
sample_data(NSQ04_RKMP_calculation_HK)$Date.Inoc <- as.factor(sample_data(NSQ04_RKMP_calculation_HK)$Date.Inoc)

#set the colours: with inoculum
SynCom_col <- c("#999999", "#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")
#set the colours: no inoculum
SynCom_col_ni <- c("#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")
#set the colours: plant
Plant_col <- c("#E1BE6A", "#40B0A6")

#Stacked bar plot
p = plot_bar(NSQ04_RKMP_calculation_HK, "Genus", facet_grid=SynCom.ID ~ Genotype, fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#Without inoculum
NSQ04_RKMP_calculation_HK_ni <- subset_samples(NSQ04_RKMP_calculation_HK, Genotype != "Inoculum")
sort(sample_sums(NSQ04_RKMP_calculation_HK_ni))

#Stacked bar plot
p = plot_bar(NSQ04_RKMP_calculation_HK_ni, "Genus", facet_grid=SynCom.ID ~ Plant, fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#alpha diversity calculation
#https://www.nature.com/articles/s41598-024-77864-y
p = plot_richness(NSQ04_RKMP_calculation_HK, x="SynCom.ID", color="Genotype", measures=c("Observed", "Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#remove Observed, this is just a confirmation
p = plot_richness(NSQ04_RKMP_calculation_HK, x="SynCom.ID", color="Genotype", measures=c("Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#without inoculum
p = plot_richness(NSQ04_RKMP_calculation_HK_ni, x="SynCom.ID", color="Plant", measures=c("Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = Plant_col)
p

####PAPER####
#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
NSQ04_RKMP_calculation_HK.ord <- ordinate(NSQ04_RKMP_calculation_HK, "NMDS", "bray")
#plot
p = plot_ordination(NSQ04_RKMP_calculation_HK, NSQ04_RKMP_calculation_HK.ord, type="samples", color="Genotype", shape="SynCom.ID") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)+
  scale_shape_discrete(name = "Assembly Method", labels = c("By OD", "By CFU"))+ theme(legend.position="bottom")
p

#unconstrained ordination - no inoculum
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
NSQ04_RKMP_calculation_HK_ni.ord <- ordinate(NSQ04_RKMP_calculation_HK_ni, "NMDS", "bray")
#plot
p = plot_ordination(NSQ04_RKMP_calculation_HK_ni, NSQ04_RKMP_calculation_HK_ni.ord, type="samples", color="Plant", shape="SynCom.ID") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = Plant_col)
p

#Constrained for SynCom effect
NSQ04_RKMP_calculation_HK_S.cap <- ordinate(NSQ04_RKMP_calculation_HK, "CAP", "bray", ~ SynCom.ID)
#Plot
a = plot_ordination(NSQ04_RKMP_calculation_HK, NSQ04_RKMP_calculation_HK_S.cap, type="samples", color="Genotype", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Constrained for SynCom effect - no inoculum
NSQ04_RKMP_calculation_HK_ni_S.cap <- ordinate(NSQ04_RKMP_calculation_HK_ni, "CAP", "bray", ~ SynCom.ID)
#Plot
a = plot_ordination(NSQ04_RKMP_calculation_HK_ni, NSQ04_RKMP_calculation_HK_ni_S.cap, type="samples", color="Plant", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = Plant_col)
a

#Constrained for 'Genotype' effect
NSQ04_RKMP_calculation_HK_G.cap <- ordinate(NSQ04_RKMP_calculation_HK, "CAP", "bray", ~ Genotype)
#Plot
a = plot_ordination(NSQ04_RKMP_calculation_HK, NSQ04_RKMP_calculation_HK_G.cap, type="samples", color="Genotype", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Constrained for plant effect - no inoculum
NSQ04_RKMP_calculation_HK_ni_P.cap <- ordinate(NSQ04_RKMP_calculation_HK_ni, "CAP", "bray", ~ Plant)
#Plot
a = plot_ordination(NSQ04_RKMP_calculation_HK_ni, NSQ04_RKMP_calculation_HK_ni_P.cap, type="samples", color="Plant", shape="SynCom.ID") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = Plant_col)
a

#Stats - Permanova
BC_bacteria  <- phyloseq::distance(NSQ04_RKMP_calculation_HK_ni, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Plant * SynCom.ID, data= as.data.frame(as.matrix(sample_data(NSQ04_RKMP_calculation_HK_ni))), permutations = 5000)
Stat

#Deseq
#extract count data 
NSQ04_RKMP_genotype_HK_ni <- otu_table(NSQ04_RKMP_calculation_HK_ni)
countData = as.data.frame(NSQ04_RKMP_genotype_HK_ni)
colnames(NSQ04_RKMP_genotype_HK_ni)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_RKMP_calculation_HK_ni)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_RKMP_HK_ni_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)

#execute the differential count analysis with the function DESeq 
NSQ04_RKMP_HK_ni_test <- DESeq(NSQ04_RKMP_HK_ni_cds, fitType="local", betaPrior = FALSE) 

#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_HK_ni_Ba <- results(NSQ04_RKMP_HK_ni_test, contrast = c("Genotype", "Bulk", "Barke")) 

#extract  ASVs whose adjusted p.value in a given comparison is below 0.05 [if any]
rhizo_genotype_RKMP_HK_ni_Ba_FDR005 <- rhizo_genotype_RKMP_HK_ni_Ba[(rownames(rhizo_genotype_RKMP_HK_ni_Ba)[which(rhizo_genotype_RKMP_HK_ni_Ba$padj <0.05)]), ]

#inspect the generated file
rhizo_genotype_RKMP_HK_ni_Ba_FDR005

#Repeat for old only

#Subset for old
NSQ04_RKMP_HK_Old <- subset_samples(NSQ04_RKMP_calculation_HK_ni, SynCom.ID == "RKMP.Old")
sort(sample_sums(NSQ04_RKMP_HK_Old))

#Repeat for New only

#Subset for New
NSQ04_RKMP_HK_New <- subset_samples(NSQ04_RKMP_calculation_HK_ni, SynCom.ID == "RKMP.New")
sort(sample_sums(NSQ04_RKMP_HK_New))
#extract count data 
NSQ04_RKMP_genotype_HK_New <- otu_table(NSQ04_RKMP_HK_New)
countData = as.data.frame(NSQ04_RKMP_genotype_HK_New)
colnames(NSQ04_RKMP_genotype_HK_New)
#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_RKMP_HK_New)))
rownames(colData)
#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_RKMP_HK_New_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)
#execute the differential count analysis with the function DESeq 
NSQ04_RKMP_HK_New_test <- DESeq(NSQ04_RKMP_HK_New_cds, fitType="local", betaPrior = FALSE) 
#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_HK_New_Ba <- results(NSQ04_RKMP_HK_New_test, contrast = c("Genotype", "Bulk", "Barke")) 
#extract  ASVs whose adjusted p.value in a given comparison is below 0.05 [if any]
rhizo_genotype_RKMP_HK_New_Ba_FDR005 <- rhizo_genotype_RKMP_HK_New_Ba[(rownames(rhizo_genotype_RKMP_HK_New_Ba)[which(rhizo_genotype_RKMP_HK_New_Ba$padj <0.05)]), ]
#inspect the generated file
rhizo_genotype_RKMP_HK_New_Ba_FDR005

#Try comparing old and new
#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_RKMP_HK_ni_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~SynCom.ID)

#execute the differential count analysis with the function DESeq 
NSQ04_RKMP_HK_ni_test <- DESeq(NSQ04_RKMP_HK_ni_cds, fitType="local", betaPrior = FALSE) 

#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_HK_ni_Ba <- results(NSQ04_RKMP_HK_ni_test, contrast = c("SynCom.ID", "RKMP.Old", "RKMP.New")) 

#extract  ASVs whose adjusted p.value in a given comparison is below 0.05 [if any]
rhizo_genotype_RKMP_HK_ni_Ba_FDR005 <- rhizo_genotype_RKMP_HK_ni_Ba[(rownames(rhizo_genotype_RKMP_HK_ni_Ba)[which(rhizo_genotype_RKMP_HK_ni_Ba$padj <0.05)]), ]

#inspect the generated file
rhizo_genotype_RKMP_HK_ni_Ba_FDR005

#Subset for Planted
NSQ04_RKMP_HK_New_P <- subset_samples(NSQ04_RKMP_calculation_HK_ni, Plant == "Planted")
sort(sample_sums(NSQ04_RKMP_HK_New_P))
#extract count data 
NSQ04_RKMP_genotype_HK_New_P <- otu_table(NSQ04_RKMP_HK_New_P)
countData = as.data.frame(NSQ04_RKMP_genotype_HK_New_P)
colnames(NSQ04_RKMP_genotype_HK_New_P)
#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_RKMP_HK_New_P)))
rownames(colData)
#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_RKMP_HK_New_P_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~SynCom.ID)
#execute the differential count analysis with the function DESeq 
NSQ04_RKMP_HK_New_P_test <- DESeq(NSQ04_RKMP_HK_New_P_cds, fitType="local", betaPrior = FALSE) 
#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_RKMP_HK_New_P_Ba <- results(NSQ04_RKMP_HK_New_P_test, contrast = c("SynCom.ID", "RKMP.Old", "RKMP.New")) 
#extract  ASVs whose adjusted p.value in a given comparison is below 0.05 [if any]
rhizo_genotype_RKMP_HK_New_P_Ba_FDR005 <- rhizo_genotype_RKMP_HK_New_P_Ba[(rownames(rhizo_genotype_RKMP_HK_New_P_Ba)[which(rhizo_genotype_RKMP_HK_New_P_Ba$padj <0.05)]), ]
#inspect the generated file
rhizo_genotype_RKMP_HK_New_P_Ba_FDR005

#########################################################################
# work with PBS samples and set criteria to work with unevenly sequenced samples
# a) keep ASV with more than 10 reads in at least 20% of the samples
# b) rarefy at even sequencing depth
# c) agglomerate at genus level
#########################################################################

#Inspect data
sort(sample_sums(NSQ04_RKMP_PBS))

#set an arbitrary threshold of 10 reads in 20% of samples
NSQ04_RKMP_PBS_threshold = filter_taxa(NSQ04_RKMP_PBS, function(x) sum(x > 10) > (0.2 *length(x)), TRUE)
NSQ04_RKMP_PBS_threshold

#inspect read distribution upon thresholding
sort(sample_sums(NSQ04_RKMP_PBS_threshold))

#agglomerate at genus level
NSQ04_RKMP_PBS_genus <- tax_glom(NSQ04_RKMP_PBS_threshold, taxrank="Genus")

#inspect the files
NSQ04_RKMP_PBS_genus
sort(sample_sums(NSQ04_RKMP_PBS_genus))

#rarefy at 260
NSQ04_RKMP_PBS_threshold_180 <- rarefy_even_depth(NSQ04_RKMP_PBS_genus, 260)

#data visualisation
p = plot_bar(NSQ04_RKMP_PBS_threshold_180, "Genus", facet_grid=SynCom.ID ~ Genotype, fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#save the file for the reproducibility of the code
#saveRDS(NSQ04_RKMP_PBS_threshold_180, file = "NSQ04_RKMP_PBS_calculation.rds")
NSQ04_RKMP_calculation_PBS <-readRDS("NSQ04_RKMP_PBS_calculation.rds")
NSQ04_RKMP_calculation_PBS

#####
#PBS#
#####

#Stacked bar plot
p = plot_bar(NSQ04_RKMP_calculation_PBS, "Genus", fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#set the colours: just Barke
SynCom_barke <- c("#0072B2")

#alpha diversity calculation
#https://www.nature.com/articles/s41598-024-77864-y
p = plot_richness(NSQ04_RKMP_calculation_PBS, x="SynCom.ID", color="Genotype", measures=c("Observed", "Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_barke)
p

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
NSQ04_RKMP_calculation_PBS.ord <- ordinate(NSQ04_RKMP_calculation_PBS, "NMDS", "bray")
#plot
p = plot_ordination(NSQ04_RKMP_calculation_PBS, NSQ04_RKMP_calculation_PBS.ord, type="samples", color="Genotype") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_barke)
p
