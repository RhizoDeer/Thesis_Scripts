#################Katie thesis 33 KOs
#Hybrid assembly only
KO_hybrid <- KO_comparison_map[KO_comparison_map$Assembly == "hybrid_26", ]
KO_counts_hybrid_t <- as.data.frame(t(KO_counts_hybrid))
row.names(KO_counts_hybrid_t) <- colnames(KO_counts_hybrid)
colnames(KO_counts_hybrid_t) <- row.names(KO_counts_hybrid)

#create the phyloseq object
counts = otu_table(KO_counts_hybrid_t, taxa_are_rows = TRUE)
samples = sample_data(KO_hybrid)

KO_SynCom_data <- phyloseq(counts, samples)
KO_SynCom_data

#list of 33 KOs
#create a phyloseq object with Pseudomonas, Stenotrophomonas, Pedobacter
KO_Pseudomonas <- subset_samples(KO_SynCom_data, Genus == "Pseudomonas")
KO_Stenotrophomonas <- subset_samples(KO_SynCom_data, Genus == "Stenotrophomonas")
KO_Pedobacter <- subset_samples(KO_SynCom_data, Genus == "Pedobacter")
#remove KOs wiht 0 counts
KO_Pseudomonas <- prune_taxa(taxa_sums(KO_Pseudomonas) > 0, KO_Pseudomonas)
KO_Stenotrophomonas <- prune_taxa(taxa_sums(KO_Stenotrophomonas) > 0, KO_Stenotrophomonas)
KO_Pedobacter <- prune_taxa(taxa_sums(KO_Pedobacter) > 0, KO_Pedobacter)
#identify the 33 conserved KOs - step1 intersection of the three genera
KO_plant_1 <- intersect(taxa_names(KO_Pseudomonas), taxa_names(KO_Stenotrophomonas))
KO_plant_1 <- intersect(KO_plant_1, taxa_names(KO_Pedobacter))
length(KO_plant_1)

#identify the 33 conserved KOs - step2 remove the background KOs [conserved between the three genera and the other "non plant" ones]
#Identify the no plants KOs
KO_no_plant <- subset_samples(KO_SynCom_data, Phylum != "Proteobacteria")
KO_no_plant <- subset_samples(KO_no_plant, Genus != "Pedobacter")
#remove KOs wiht 0 counts
KO_no_plant  <- prune_taxa(taxa_sums(KO_no_plant ) > 0, KO_no_plant)
KO_no_plant
#uniquely conserved by the three genera
Plant_KOs <- setdiff(KO_plant_1, taxa_names(KO_no_plant))
 
#Actual KOs [copied Plant entries in Plant_KOs here https://www.kegg.jp/kegg/ko.html, the "enter the KO number allows for multiple entries]
#K00209  fabV, ter; enoyl-[acyl-carrier protein] reductase / trans-2-enoyl-CoA reductase (NAD+) [EC:1.3.1.9 1.3.1.44]
#K00228  CPOX, hemF; coproporphyrinogen III oxidase [EC:1.3.3.3]
#K00500  phhA, PAH; phenylalanine-4-hydroxylase [EC:1.14.16.1]
#K00507  SCD, desC; stearoyl-CoA desaturase (Delta-9 desaturase) [EC:1.14.19.1]
#K01070  frmB, ESD, fghA; S-formylglutathione hydrolase [EC:3.1.2.12]
#K01241  amn; AMP nucleosidase [EC:3.2.2.4]
#K01729  algL; poly(beta-D-mannuronate) lyase [EC:4.2.2.3]
#K02193  ccmA; heme exporter protein A [EC:7.6.2.5]
#K02194  ccmB; heme exporter protein B
#K02195  ccmC; heme exporter protein C
#K02197  ccmE; cytochrome c-type biogenesis protein CcmE
#K02198  ccmF; cytochrome c-type biogenesis protein CcmF
#K02454  gspE; general secretion pathway protein E [EC:7.4.2.8]
#K02653  pilC; type IV pilus assembly protein PilC
#K02655  pilE; type IV pilus assembly protein PilE
#K02666  pilQ; type IV pilus assembly protein PilQ
#K03592  pmbA; PmbA protein
#K03778  ldhA; D-lactate dehydrogenase [EC:1.1.1.28]
#K04046  yegD; hypothetical chaperone protein
#K05367  pbpC; penicillin-binding protein 1C [EC:2.4.99.28]
#K05802  mscK, kefA, aefA; potassium-dependent mechanosensitive channel
#K06179  rluC; 23S rRNA pseudouridine955/2504/2580 synthase [EC:5.4.99.24]
#K06894  yfhM; alpha-2-macroglobulin
#K07506  K07506; AraC family transcriptional regulator
#K07679  evgS, bvgS; two-component system, NarL family, sensor histidine kinase EvgS [EC:2.7.13.3]
#K08482  kaiC; circadian clock protein KaiC
#K09915  K09915; uncharacterized protein
#K09933  mtfA; MtfA peptidase
#K09950  K09950; uncharacterized protein
#K11719  lptC; lipopolysaccharide export system protein LptC
#K16092  btuB; vitamin B12 transporter
#K18139  oprM, emhC, ttgC, cusC, adeK, smeF, mtrE, cmeC, gesC; outer membrane protein, multidrug efflux system
#K18300  oprN; outer membrane protein, multidrug efflux system