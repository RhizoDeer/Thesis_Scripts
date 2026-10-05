#####################################################################################
#Figure KO comparison old vs new bacterial genome sequencing
####################################################################################
#############################################################
#
# Ref to the ARTICLE 
# 
#  
#  Revision 08/26
#  c.arnton@dundee.ac.uk 
#  d.bulgarelli@dundee.ac.uk 
#
#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################
library("phyloseq")
library("vegan")
library("ggplot2")
library("UpSetR")
library("viridis")
#############################################################
#set working directory: DB cpu
setwd("/cluster/db/R_shared/MicrobesNG_0426/")
#set working directory: Katie home
setwd("C:/Users/catar/Documents/Dundee PhD/Thesis/Chapter 3 Analysis")
getwd()

#import the constituent components
#KO counts short reads
KO_counts_short <- read.delim("ko_SynCom_strains_Short.txt", row.names =1)
dim(KO_counts_short)

#KO counts hybrid
KO_counts_hybrid <- read.delim("ko_SynCom_strains_Hybrid.txt", row.names =1)
dim(KO_counts_hybrid)

#proportion of conserved KOs
ratio <- length(intersect(colnames(KO_counts_short),colnames(KO_counts_hybrid)))/length(unique(union(colnames(KO_counts_short),colnames(KO_counts_hybrid)))) * 100
ratio 

#samples
KO_comparison_map <- read.delim("ko_comparison_mapping_2.txt", row.names =1)

#transpose the datasets
KO_counts_short_t <- as.data.frame(t(KO_counts_short))
row.names(KO_counts_short_t) <- colnames(KO_counts_short)
colnames(KO_counts_short_t) <- row.names(KO_counts_short)
KO_counts_hybrid_t <- as.data.frame(t(KO_counts_hybrid))
row.names(KO_counts_hybrid_t) <- colnames(KO_counts_hybrid)
colnames(KO_counts_hybrid_t) <- row.names(KO_counts_hybrid)

#merge the ko counts table: note the have unequal columns
KO_list <- unique(c(row.names(KO_counts_short_t), (row.names(KO_counts_hybrid_t))))
length(KO_list)

#different method for conserved KOs
length(intersect(colnames(KO_counts_short),colnames(KO_counts_hybrid)))/length(KO_list)*100

#configure the data frame
KO_counts_short_merging <-  as.data.frame(KO_counts_short_t[KO_list, ])
KO_counts_hybrid_merging <-  as.data.frame(KO_counts_hybrid_t[KO_list, ])
#merge the dataset as now they have the same rows
KO_comparison_counts <-  cbind(KO_counts_short_merging, KO_counts_hybrid_merging)
#set to 0 rows with NA, i.e., the one missing in one dataset
KO_comparison_counts[is.na(KO_comparison_counts)] <- 0
dim(KO_comparison_counts)

#create the phyloseq object
counts = otu_table(KO_comparison_counts, taxa_are_rows = TRUE)
samples = sample_data(KO_comparison_map)

KO_SynCom_data <- phyloseq(counts, samples)
KO_SynCom_data

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
KO_SynCom_data.ord <- ordinate(KO_SynCom_data, "NMDS", "bray")
#plot
p = plot_ordination(KO_SynCom_data, KO_SynCom_data.ord, type="samples", color="Genus", shape = "Assembly") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_viridis(discrete=TRUE)
p

####Thesis####
#constrain for both assembly and genus
KO_SynCom.cap <- ordinate(KO_SynCom_data, "CAP", "bray", ~ Assembly * Genus)
#ggplots function to increase effectivness of the visualisation
p = plot_ordination(KO_SynCom_data, KO_SynCom.cap,  color="Genus", shape="Assembly") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_viridis(discrete=TRUE)
p = p + scale_shape_discrete(name = "Assembly Method", labels = c("Hybrid", "Illumina"))
p + theme(legend.position="right") 
#all near their partner with the other assembly, stenotrophomonas slightly further away than rest

#Just genus
KO_SynCom.cap <- ordinate(KO_SynCom_data, "CAP", "bray", ~ Genus)
p = plot_ordination(KO_SynCom_data, KO_SynCom.cap,  color="Genus", shape="Assembly") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_viridis(discrete=TRUE)
p = p + scale_shape_discrete(name = "Assembly Method", labels = c("Hybrid", "Illumina"))
p + theme(legend.position="right") 

#this calculation should use the same distance used to build the graphical output
BC_KO  <- phyloseq::distance(KO_SynCom_data, "bray")
BC_KO

#here we can use a formula ANOVA-like to identify factors of interest
Stat_KO <- adonis2(BC_KO ~ Assembly * Genus, data= as.data.frame(as.matrix(sample_data(KO_SynCom_data))), permutations = 5000)
Stat_KO

#move to phylum level
Stat_KO <- adonis2(BC_KO ~ Assembly * Phylum, data= as.data.frame(as.matrix(sample_data(KO_SynCom_data))), permutations = 5000)
Stat_KO

#run the same analysis with the gram variable
Stat_KO <- adonis2(BC_KO ~ Assembly * Gram, data= as.data.frame(as.matrix(sample_data(KO_SynCom_data))), permutations = 5000)
Stat_KO

#identify KO differentially coded - split the database for phyla and run an upsetR plot
KO_SynCom_data_gram_plus <- subset_samples(KO_SynCom_data, Gram == "positive")
KO_SynCom_data_gram_minus <- subset_samples(KO_SynCom_data, Gram == "negative")

#remove rows with 0 counts
KO_SynCom_data_gram_plus <- prune_taxa(taxa_sums(KO_SynCom_data_gram_plus) > 0, KO_SynCom_data_gram_plus)
KO_SynCom_data_gram_minus <- prune_taxa(taxa_sums(KO_SynCom_data_gram_minus) > 0, KO_SynCom_data_gram_minus)

#gram negative enriched
gram_miuns_function <- setdiff(as.vector(taxa_names(KO_SynCom_data_gram_minus)), as.vector(taxa_names(KO_SynCom_data_gram_plus)))
gram_plus_function <- setdiff(as.vector(taxa_names(KO_SynCom_data_gram_plus)), as.vector(taxa_names(KO_SynCom_data_gram_minus)))

#export the files
write(gram_miuns_function, file = "KO_minus_function.txt")
write(gram_plus_function, file = "KO_plus_function.txt")

#plotting KOs differentially coded among taxa => collapse at phylum level
KO_SynCom_data_Actinobacteria <- subset_samples(KO_SynCom_data, Phylum == "Actinobacteria")
KO_SynCom_data_Actinobacteria <- prune_taxa(taxa_sums(KO_SynCom_data_Actinobacteria) > 0, KO_SynCom_data_Actinobacteria)
KO_SynCom_data_Bacteroidetes <- subset_samples(KO_SynCom_data, Phylum == "Bacteroidetes")
KO_SynCom_data_Bacteroidetes <- prune_taxa(taxa_sums(KO_SynCom_data_Bacteroidetes) > 0, KO_SynCom_data_Bacteroidetes)
KO_SynCom_data_Firmicutes <- subset_samples(KO_SynCom_data, Phylum == "Firmicutes")
KO_SynCom_data_Firmicutes <- prune_taxa(taxa_sums(KO_SynCom_data_Firmicutes) > 0, KO_SynCom_data_Firmicutes)
KO_SynCom_data_Proteobacteria <- subset_samples(KO_SynCom_data, Phylum == "Proteobacteria")
KO_SynCom_data_Proteobacteria <- prune_taxa(taxa_sums(KO_SynCom_data_Proteobacteria) > 0, KO_SynCom_data_Proteobacteria)

#create a new dataset with unique KOs per Phylum
#Actinobacteria
Actinobacteria_unique_counts <- as.data.frame(taxa_sums(KO_SynCom_data_Actinobacteria))
dim(Actinobacteria_unique_counts)
colnames(Actinobacteria_unique_counts) <- c("KOs_Actinobacteria")
Actinobacteria_unique_counts[Actinobacteria_unique_counts > 1] <- 1
dim(Actinobacteria_unique_counts)
#Bacteroidetes
Bacteroidetes_unique_counts <- as.data.frame(taxa_sums(KO_SynCom_data_Bacteroidetes))
dim(Bacteroidetes_unique_counts)
colnames(Bacteroidetes_unique_counts) <- c("KOs_Bacteroidetes")
Bacteroidetes_unique_counts[Bacteroidetes_unique_counts > 1] <- 1
dim(Bacteroidetes_unique_counts)
#Firmicutes
Firmicutes_unique_counts <- as.data.frame(taxa_sums(KO_SynCom_data_Firmicutes))
dim(Firmicutes_unique_counts)
colnames(Firmicutes_unique_counts) <- c("KOs_Firmicutes")
Firmicutes_unique_counts[Firmicutes_unique_counts > 1] <- 1
dim(Firmicutes_unique_counts)
#Proteobacteria
Proteobacteria_unique_counts <- as.data.frame(taxa_sums(KO_SynCom_data_Proteobacteria))
dim(Proteobacteria_unique_counts)
colnames(Proteobacteria_unique_counts) <- c("KOs_Proteobacteria")
Proteobacteria_unique_counts[Proteobacteria_unique_counts > 1] <- 1
dim(Proteobacteria_unique_counts)

#merging the dataset
#Prior combining the dataset, we need to account for ASV unevenly distributed (i.e., enriched in one compartment not in others)
#Actinobacteria
Actinobacteria_unique_merging <- as.data.frame(Actinobacteria_unique_counts[KO_list, ])
colnames(Actinobacteria_unique_merging) <- c("KOs_Actinobacteria")
row.names(Actinobacteria_unique_merging) <- as.vector(KO_list)
#Firmicutes
Firmicutes_unique_merging <- as.data.frame(Firmicutes_unique_counts[KO_list, ])
colnames(Firmicutes_unique_merging) <- c("KOs_Firmicutes")
row.names(Firmicutes_unique_merging) <- as.vector(KO_list)
#Bacteroidetes
Bacteroidetes_unique_merging <- as.data.frame(Bacteroidetes_unique_counts[KO_list, ])
colnames(Bacteroidetes_unique_merging) <- c("KOs_Bacteroidetes")
row.names(Bacteroidetes_unique_merging) <- as.vector(KO_list)
#Proteobacteria
Proteobacteria_unique_merging <- as.data.frame(Proteobacteria_unique_counts[KO_list, ])
colnames(Proteobacteria_unique_merging) <- c("KOs_Proteobacteria")
row.names(Proteobacteria_unique_merging) <- as.vector(KO_list)

#Merge the datasets
Strains_KOs <- cbind(Actinobacteria_unique_merging, Bacteroidetes_unique_merging)
Strains_KOs <- cbind(Strains_KOs, Firmicutes_unique_merging)
Strains_KOs <- cbind(Strains_KOs, Proteobacteria_unique_merging)
#set NA to 0: as some ASVs won't be present in certain microhabitats and NA is a non numerical character
Strains_KOs[is.na(Strains_KOs)] <- 0
dim(Strains_KOs)

#Plot
plot_UPSET_KOs_phyla <-upset(Strains_KOs, sets = c("KOs_Actinobacteria", "KOs_Bacteroidetes", "KOs_Firmicutes", "KOs_Proteobacteria"), sets.bar.color = "#56B4E9",
                               order.by = "freq", sets.x.label = "KOs coded", mainbar.y.label = "Intersection",)
plot_UPSET_KOs_phyla

#An upsetR plot using two completely independent assemblies makes no biological sense, repeat with just hybrid

#Make phyloseq object
KO_comparison_map_h <- read.delim("ko_comparison_mapping_hybrid.txt", row.names =1)
counts = otu_table(KO_counts_hybrid_t, taxa_are_rows = TRUE)
samples = sample_data(KO_comparison_map_h)
KO_SynCom_data_h <- phyloseq(counts, samples)
KO_SynCom_data_h

#unconstrained ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
KO_SynCom_data_h.ord <- ordinate(KO_SynCom_data_h, "NMDS", "bray")
#plot
p = plot_ordination(KO_SynCom_data_h, KO_SynCom_data_h.ord, type="samples", color="Genus") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_viridis(discrete=TRUE)
p

#constrain for genus
KO_SynCom.cap <- ordinate(KO_SynCom_data_h, "CAP", "bray", ~Genus)
#ggplots function to increase effectivness of the visualisation
p = plot_ordination(KO_SynCom_data_h, KO_SynCom.cap,  color="Genus") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_viridis(discrete=TRUE)
p + theme(legend.position="right") 

#plotting KOs differentially coded among taxa => collapse at genus level
KO_SynCom_data_h_Arthrobacter <- subset_samples(KO_SynCom_data_h, Genus == "Arthrobacter")
KO_SynCom_data_h_Arthrobacter <- prune_taxa(taxa_sums(KO_SynCom_data_h_Arthrobacter) > 0, KO_SynCom_data_h_Arthrobacter)
KO_SynCom_data_h_Chryseobacterium <- subset_samples(KO_SynCom_data_h, Genus == "Chryseobacterium")
KO_SynCom_data_h_Chryseobacterium <- prune_taxa(taxa_sums(KO_SynCom_data_h_Chryseobacterium) > 0, KO_SynCom_data_h_Chryseobacterium)
KO_SynCom_data_h_Pedobacter <- subset_samples(KO_SynCom_data_h, Genus == "Pedobacter")
KO_SynCom_data_h_Pedobacter <- prune_taxa(taxa_sums(KO_SynCom_data_h_Pedobacter) > 0, KO_SynCom_data_h_Pedobacter)
KO_SynCom_data_h_Peribacillus <- subset_samples(KO_SynCom_data_h, Genus == "Peribacillus")
KO_SynCom_data_h_Peribacillus <- prune_taxa(taxa_sums(KO_SynCom_data_h_Peribacillus) > 0, KO_SynCom_data_h_Peribacillus)
KO_SynCom_data_h_Priestia <- subset_samples(KO_SynCom_data_h, Genus == "Priestia")
KO_SynCom_data_h_Priestia <- prune_taxa(taxa_sums(KO_SynCom_data_h_Priestia) > 0, KO_SynCom_data_h_Priestia)
KO_SynCom_data_h_Pseudomonas <- subset_samples(KO_SynCom_data_h, Genus == "Pseudomonas")
KO_SynCom_data_h_Pseudomonas <- prune_taxa(taxa_sums(KO_SynCom_data_h_Pseudomonas) > 0, KO_SynCom_data_h_Pseudomonas)
KO_SynCom_data_h_Rhodococcus <- subset_samples(KO_SynCom_data_h, Genus == "Rhodococcus")
KO_SynCom_data_h_Rhodococcus <- prune_taxa(taxa_sums(KO_SynCom_data_h_Rhodococcus) > 0, KO_SynCom_data_h_Rhodococcus)
KO_SynCom_data_h_Stenotrophomonas <- subset_samples(KO_SynCom_data_h, Genus == "Stenotrophomonas")
KO_SynCom_data_h_Stenotrophomonas <- prune_taxa(taxa_sums(KO_SynCom_data_h_Stenotrophomonas) > 0, KO_SynCom_data_h_Stenotrophomonas)

#create a new dataset with unique KOs per genus
#Arthrobacter
Arthrobacter_unique_counts <- as.data.frame(taxa_sums(KO_SynCom_data_h_Arthrobacter))
dim(Arthrobacter_unique_counts)
colnames(Arthrobacter_unique_counts) <- c("KOs_Arthrobacter")
Arthrobacter_unique_counts[Arthrobacter_unique_counts > 1] <- 1
dim(Arthrobacter_unique_counts)
#Chryseobacterium
Chryseobacterium_unique_counts <- as.data.frame(taxa_sums(KO_SynCom_data_h_Chryseobacterium))
dim(Chryseobacterium_unique_counts)
colnames(Chryseobacterium_unique_counts) <- c("KOs_Chryseobacterium")
Chryseobacterium_unique_counts[Chryseobacterium_unique_counts > 1] <- 1
dim(Chryseobacterium_unique_counts)
#Pedobacter
Pedobacter_unique_counts <- as.data.frame(taxa_sums(KO_SynCom_data_h_Pedobacter))
dim(Pedobacter_unique_counts)
colnames(Pedobacter_unique_counts) <- c("KOs_Pedobacter")
Pedobacter_unique_counts[Pedobacter_unique_counts > 1] <- 1
dim(Pedobacter_unique_counts)
#Peribacillus
Peribacillus_unique_counts <- as.data.frame(taxa_sums(KO_SynCom_data_h_Peribacillus))
dim(Peribacillus_unique_counts)
colnames(Peribacillus_unique_counts) <- c("KOs_Peribacillus")
Peribacillus_unique_counts[Peribacillus_unique_counts > 1] <- 1
dim(Peribacillus_unique_counts)
#Priestia
Priestia_unique_counts <- as.data.frame(taxa_sums(KO_SynCom_data_h_Priestia))
dim(Priestia_unique_counts)
colnames(Priestia_unique_counts) <- c("KOs_Priestia")
Priestia_unique_counts[Priestia_unique_counts > 1] <- 1
dim(Priestia_unique_counts)
#Pseudomonas
Pseudomonas_unique_counts <- as.data.frame(taxa_sums(KO_SynCom_data_h_Pseudomonas))
dim(Pseudomonas_unique_counts)
colnames(Pseudomonas_unique_counts) <- c("KOs_Pseudomonas")
Pseudomonas_unique_counts[Pseudomonas_unique_counts > 1] <- 1
dim(Pseudomonas_unique_counts)
#Rhodococcus
Rhodococcus_unique_counts <- as.data.frame(taxa_sums(KO_SynCom_data_h_Rhodococcus))
dim(Rhodococcus_unique_counts)
colnames(Rhodococcus_unique_counts) <- c("KOs_Rhodococcus")
Rhodococcus_unique_counts[Rhodococcus_unique_counts > 1] <- 1
dim(Rhodococcus_unique_counts)
#Stenotrophomonas
Stenotrophomonas_unique_counts <- as.data.frame(taxa_sums(KO_SynCom_data_h_Stenotrophomonas))
dim(Stenotrophomonas_unique_counts)
colnames(Stenotrophomonas_unique_counts) <- c("KOs_Stenotrophomonas")
Stenotrophomonas_unique_counts[Stenotrophomonas_unique_counts > 1] <- 1
dim(Stenotrophomonas_unique_counts)

#merging the dataset
#Prior combining the dataset, we need to account for ASV unevenly distributed (i.e., enriched in one compartment not in others)
KO_list_h <- row.names(KO_counts_hybrid_t)
length(KO_list_h)
#Arthrobacter
Arthrobacter_unique_merging <- as.data.frame(Arthrobacter_unique_counts[KO_list_h, ])
colnames(Arthrobacter_unique_merging) <- c("KOs_Arthrobacter")
row.names(Arthrobacter_unique_merging) <- as.vector(KO_list_h)
#Chryseobacterium
Chryseobacterium_unique_merging <- as.data.frame(Chryseobacterium_unique_counts[KO_list_h, ])
colnames(Chryseobacterium_unique_merging) <- c("KOs_Chryseobacterium")
row.names(Chryseobacterium_unique_merging) <- as.vector(KO_list_h)
#Pedobacter
Pedobacter_unique_merging <- as.data.frame(Pedobacter_unique_counts[KO_list_h, ])
colnames(Pedobacter_unique_merging) <- c("KOs_Pedobacter")
row.names(Pedobacter_unique_merging) <- as.vector(KO_list_h)
#Peribacillus
Peribacillus_unique_merging <- as.data.frame(Peribacillus_unique_counts[KO_list_h, ])
colnames(Peribacillus_unique_merging) <- c("KOs_Peribacillus")
row.names(Peribacillus_unique_merging) <- as.vector(KO_list_h)
#Priestia
Priestia_unique_merging <- as.data.frame(Priestia_unique_counts[KO_list_h, ])
colnames(Priestia_unique_merging) <- c("KOs_Priestia")
row.names(Priestia_unique_merging) <- as.vector(KO_list_h)
#Pseudomonas
Pseudomonas_unique_merging <- as.data.frame(Pseudomonas_unique_counts[KO_list_h, ])
colnames(Pseudomonas_unique_merging) <- c("KOs_Pseudomonas")
row.names(Pseudomonas_unique_merging) <- as.vector(KO_list_h)
#Rhodococcus
Rhodococcus_unique_merging <- as.data.frame(Rhodococcus_unique_counts[KO_list_h, ])
colnames(Rhodococcus_unique_merging) <- c("KOs_Rhodococcus")
row.names(Rhodococcus_unique_merging) <- as.vector(KO_list_h)
#Stenotrophomonas
Stenotrophomonas_unique_merging <- as.data.frame(Stenotrophomonas_unique_counts[KO_list_h, ])
colnames(Stenotrophomonas_unique_merging) <- c("KOs_Stenotrophomonas")
row.names(Stenotrophomonas_unique_merging) <- as.vector(KO_list_h)

#Merge the datasets
Strains_KOs <- cbind(Arthrobacter_unique_merging, Chryseobacterium_unique_merging)
Strains_KOs <- cbind(Strains_KOs, Pedobacter_unique_merging)
Strains_KOs <- cbind(Strains_KOs, Peribacillus_unique_merging)
Strains_KOs <- cbind(Strains_KOs, Priestia_unique_merging)
Strains_KOs <- cbind(Strains_KOs, Pseudomonas_unique_merging)
Strains_KOs <- cbind(Strains_KOs, Rhodococcus_unique_merging)
Strains_KOs <- cbind(Strains_KOs, Stenotrophomonas_unique_merging)
#set NA to 0: as some ASVs won't be present in certain microhabitats and NA is a non numerical character
Strains_KOs[is.na(Strains_KOs)] <- 0
dim(Strains_KOs)

#Plot
plot_UPSET_KOs_phyla <-upset(Strains_KOs, sets = c("KOs_Arthrobacter", "KOs_Chryseobacterium", "KOs_Pedobacter", "KOs_Peribacillus", "KOs_Priestia", "KOs_Pseudomonas", "KOs_Rhodococcus", "KOs_Stenotrophomonas"), sets.bar.color = "#56B4E9",
                             order.by = "freq", sets.x.label = "KOs coded", mainbar.y.label = "Intersection",)
plot_UPSET_KOs_phyla
