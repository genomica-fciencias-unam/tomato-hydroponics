# Load packages and phyloseq data


```R
# Loading packages
library(phyloseq)
library(vegan)
library(ggplot2)
library(VennDiagram)
library(magrittr)
library(ape)
library(ggsignif)
library(DESeq2)
library(limma)
library(tidyverse)
library(scales)
library(reshape2)
library(dplyr)
library(UpSetR)
library(tibble)
```


```R
#Loading phyloseq object
otu <- as.matrix(read.table("all_proteins_global_ss_counts.tsv", header=T, row.names=1)) 
OTU = otu_table(otu, taxa_are_rows=T)
taximat = as.matrix(read.table("all_proteins_global_ss_onto.tsv", header=T, row.names=1)) 
taxi=tax_table(taximat)
tom_all= phyloseq(OTU, taxi)
sample_names(tom_all)
data= read.table("metada.csv", header=T, row.names=1, sep=",")
head(data)
sampledata = sample_data(data.frame(sample=data$sample,subid=data$subsam, type=data$type,
                                    ID=data$ID, substratum=data$substratum, 
                                    row.names=sample_names(tom_all)))
tom_all= phyloseq (OTU, sampledata, taxi)
tom_all
```

# Alpha diversity indices


```R
#plot diversity indices
ind_alfa<-plot_richness(tom_all, x="substratum", color="substratum", 
                        measures=c("Observed", "Shannon", "Simpson")) + geom_boxplot() + geom_point(size=4) +
theme_light()+ scale_colour_manual(values = c("steelblue","rosybrown4")) + 
geom_signif(test = 'wilcox.test', map_signif_level = TRUE, comparisons = list(c("hydroponic", "soil")), textsize=6) + 
theme(axis.text.x=element_text(angle=90, hjust=1))
div_table<-estimate_richness(tom_all , measures=c("Observed", "Chao1", "Shannon", "Simpson"))
ind_alfa
```

# Ordenamiento NMDS


```R
#calculate distance matrix
gp.ord <- ordinate(tom_all, "NMDS", "bray")
gp.ord
#Plot NMDS
plot_ordination(tom_all, gp.ord, color="substratum") + theme_bw () + theme(text = element_text(size = 16)) + 
geom_point(size = 4) + 
geom_text(mapping = aes(label = subid), size = 3)  + 
scale_color_manual(values = c("steelblue","rosybrown4")) +theme_light()
```

# Identify overrepresented proteins


```R
#Subset hydroponics
p_hydro<-filter_taxa(subset_samples(tom_all, substratum=="hydroponic"), function (x) {sum(x > 0) > 0}, prune=TRUE)
p_hydro
p_hydross<-filter_taxa(subset_samples(p_hydro, typy=="ss"), function (x) {sum(x > 0) > 0}, prune=TRUE)
p_hydross

#Subset soil
p_soil<-filter_taxa(subset_samples(tom_all, substratum=="soil"), function (x) {sum(x > 0) > 0}, prune=TRUE)
p_soil
p_soilss<-filter_taxa(subset_samples(p_soil, typy=="ss"), function (x) {sum(x > 0) > 0}, prune=TRUE)
p_soilss

#Save taxonomy tables
write.table((tax_table(p_hydross)), "taxt_allp_hydro_ss.tsv", sep="\t")
write.table((tax_table(p_soilss)), "taxt_allp_soil_ss.tsv", sep="\t")

#Load taxonomy tables
c_hydro <- read.table("taxt_allp_hydro_ss.tsv",header=T, row.names=NULL)
c_hydro <-as.vector(c_hydro$row.names)

c_soil <- read.table("taxt_allp_soil_ss.tsv",header=T, row.names=NULL)
c_soil <-as.vector(c_soil$row.names)

#Create a list
l_hydro_soil<-list(c_hydro, c_soil)

#Creating functions to obtain the lists for each subset
Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

#Obtain the intersection of proteins 
l_hydro_soil_s<-Intersect(l_hydro_soil)
l_hydro_soil_d<-as.data.frame(l_hydro_soil_s,row.names = NULL)

#obtener los OTUs compartidos por grown hydro and soil
subset_hydroysoil<- subset_taxa(tom_all, rownames(tax_table(tom_all)) %in% l_hydro_soil_s)
subset_hydroysoil

#Enrichement analysis of controles de hydro and soil

print("1")
alpha = 0.01

merge_phyloseq(subset_samples(subset_hydroysoil, substratum=="hydroponic"), 
               subset_samples(subset_hydroysoil,substratum=="soil"))->sar.shu.phy
sar.shu.phy

print("2")
sar.shu.ds <- phyloseq_to_deseq2(sar.shu.phy, ~substratum)
print("3")
sar.shu.ds<-DESeq(sar.shu.ds, test="Wald", fitType="local")
sar.shu.ds.res <- results(sar.shu.ds, cooksCutoff = FALSE)
sar.shu.ds.res
print("4")
sigtab.sar.shu<-sar.shu.ds.res[which(sar.shu.ds.res$padj < alpha), ]
print("5")
sigtab.sar.shu<-cbind(as(sigtab.sar.shu, "data.frame"), 
                      as(tax_table(sar.shu.phy)[rownames(sigtab.sar.shu), ], "matrix"))
print("6")
sigtab.sar.shu.x=tapply(sigtab.sar.shu$log2FoldChange, sigtab.sar.shu$protein1, 
                        function(x) max(x))
                        
print("7")
sigtab.sar.shu.x=sort(sigtab.sar.shu.x, TRUE)
print("8")
sigtab.sar.shu$Genus = factor(as.character(sigtab.sar.shu$protein1), 
                              levels=names(sigtab.sar.shu.x))              
write.csv(sigtab.sar.shu, "log2_hydroysoil_01.csv")
```


```R
#Plot significant enriched proteins
sigtab.sar.shu <- read.table("log2_hydroysoil_01.csv",header=T, row.names=NULL, sep=",")

#Significance threshold
umbral <- 0.01

#Create a volcano plot with custom colors
volcano_plot <- ggplot(sigtab.sar.shu, aes(x = -log10(pvalue), y = log2FoldChange)) +
  geom_point(aes(color = ifelse(log2FoldChange < 0, "Negative", "Positive")), size = 2, alpha = 1) +
  scale_color_manual(
    values = c("Negative" = "#4e6b98ff", "Positive" = "#997973ff"),
    guide = "none"
  ) +
  theme_light() +
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5)
  )
volcano_plot
```

# COG classification of overrepresented proteins in hydroponics


```R
#Sequences corresponding to overrepresented proteins in hydroponics were annotated with COG using the EggnogMapper online tool with default parameters
```


```R
# seqprot_names_idsm5_sobreh_001.txt contains the relationship of all overrepresented proteins between the sequence name and its ID in M5nr
# ot_rsq_r.tsv contains the relative abundance of the overrepresented proteins
# sobrer_hydrovssoil_COG_001.csv contains COG annotation by sequence name


# Load tables
seq_idm5 <- read.table("seqprot_names_idsm5_sobreh_001.txt", header=FALSE, row.names=NULL, 
                       stringsAsFactors = FALSE, check.names=FALSE)
colnames(seq_idm5) <- c("protein", "idm5")

idm5_srfs <- read.table("ot_rsq_r.tsv", header=TRUE, row.names=NULL, sep="\t")
colnames(idm5_srfs)[1] <- c("idm5")

cog_seq <- read.table("sobrer_hydrovssoil_COG_001.csv", header=TRUE, row.names=NULL, sep=",")
colnames(cog_seq)[1] <- c("protein")

# Merge the tables based on the common columns
seq_idm5_srfs <- merge(x=seq_idm5, y=idm5_srfs, by="idm5")
cog_seq_idm5_srfs <- merge(x=seq_idm5_srfs, y=cog_seq, by="protein")

# Extract relevant columns for COG abundance
cog_ab <- cog_seq_idm5_srfs[,c(35,3:34)]
cog_ab_order <- melt(cog_ab)

# Group by COG categories and variables, summarizing the values
cog_group <- cog_ab_order %>%
  group_by(COG_category, variable) %>%
  summarise(value = sum(value))

# Remove rows with zero values
cog_group_no0 <- cog_group[cog_group$value != 0,]

# Order by value in descending order
cog_group_no0_order <- cog_group_no0[order(cog_group_no0$value, decreasing=TRUE),]

# Plot relative abundance
ggplot(cog_group_no0_order, aes(x=variable, y=reorder(COG_category, value), fill=value)) + 
  geom_raster() + theme_minimal() + theme(text = element_text(size=11)) +
  scale_fill_gradient(low="#f9e1e0ff", high="#5f618fff", na.value = "white", trans = "log10") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

# Taxonomic classification of overrepresented hydroponic proteins


```R
# Subset the overrepresented proteins
sobrer_sub <- read.table("id_sobrer_hydro_001.txt", header=FALSE, row.names=NULL, 
                         stringsAsFactors = FALSE, check.names=FALSE)
p1 <- sobrer_sub$V1

# Subset the taxa using p1 and perform taxonomic agglomeration
sobrer_sub_p1 <- subset_taxa(tom_all_rsq, rownames(tax_table(tom_all_rsq)) %in% p1)
sobrer_sub_p1

# Save ontology table
onto_table <- as.data.frame(tax_table(sobrer_sub_p1))
onto_table$id <- row.names(onto_table)
row.names(onto_table) <- NULL

# Load tables 
table_1 <- read.table("sobrer_hydro_001.tsv", header=FALSE, row.names=NULL,
                      check.names=FALSE)
colnames(table_1) <- c("id", "protein", "contig_a")

table_2 <- read.table("sobrer_hydro_ann.tsv", header=FALSE, row.names=NULL,
                      check.names=FALSE, sep = "\t")
colnames(table_2) <- c("contig_b", "tax", "tax_id")

# sobrer_hydro_001.tsv contains the list of identifiers for overrepresented proteins in hydroponics
# sobrer_hydro_ann.tsv contains the taxonomic assignment for each protein identifier

# Merge tables
m1 <- cbind(table_1, table_2)
m2 <- merge(onto_table, m1, by="id")

# Count all entries
freq_percent <- m2 %>% 
  count(tax) %>%  
  mutate(perc = round(n/sum(n)*100, 2)) %>% 
  arrange(desc(perc))

# Count entries with a percentage >= 0.5
freq_percent05 <- m2 %>% 
  count(tax) %>%  
  mutate(perc = round(n/sum(n)*100, 2)) %>% 
  filter(perc >= 0.5) %>% 
  arrange(desc(perc))

# Aggregate low-frequency taxa
fres_col_low <- freq_percent %>%
  mutate(tax = ifelse(perc > 0.1, tax, "low")) %>%
  group_by(tax) %>%
  summarise(perc_sum = sum(perc))

# Create plot for aggregated data
ggplot(fres_col_low, aes(y = reorder(tax, +perc_sum), x = perc_sum)) +
  geom_bar(stat = "identity") +
  labs(x = "Taxon", y = "Frequency (%)") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(panel.background = element_blank()) +
  theme(text = element_text(size=20)) 
```

# Core calculations for soil and hydroponics


```R
#Creating functions to generate lists for each subset

Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}
```


```R
#core hydroponics
#filtering hydro samples 
hydro<-filter_taxa(subset_samples(tom_all, substratum=="hydroponic"), function (x) {sum(x > 0) > 0}, prune=TRUE)
hydro
hydro<-filter_taxa(subset_samples(hydro, type=="sample"), function (x) {sum(x > 0) > 0}, prune=TRUE)
hydro

#creating upset
entradaUpset <- otu_table(hydro)
entradaUpset[entradaUpset>0] = 1
tomateUpset <- upset(as.data.frame(entradaUpset), nsets=40, nintersects=7, order.by = "degree", number.angles = 0, point.size = 2.5, line.size = 0.5, mainbar.y.label = "Protein intersections", 
                     sets.x.label = "Proteins per sample",  text.scale = c(4,4,4,4,4,4))
tomateUpset
```


```R
#Get the core hydroponics list 
#Subset samples
hydro2 <- filter_taxa(subset_samples(hydro, sample == "hydro_2"), function(x) { sum(x > 0) > 0 }, prune = TRUE)
hydro3 <- filter_taxa(subset_samples(hydro, sample == "hydro_3"), function(x) { sum(x > 0) > 0 }, prune = TRUE)
hydro4 <- filter_taxa(subset_samples(hydro, sample == "hydro_4"), function(x) { sum(x > 0) > 0 }, prune = TRUE)

# Save proteins list
hydro2_prot <- as.data.frame(tax_table(hydro2))
hydro2_prot_v <- as.vector(row.names(hydro2_prot))
hydro3_prot <- as.data.frame(tax_table(hydro3))
hydro3_prot_v <- as.vector(row.names(hydro3_prot))
hydro4_prot <- as.data.frame(tax_table(hydro4))
hydro4_prot_v <- as.vector(row.names(hydro4_prot))

# Obtain the intersection of all proteins
h2h3h4_l <- list(hydro2_prot_v, hydro3_prot_v, hydro4_prot_v)
i_h2h3h4 <- Intersect(h2h3h4_l)
i_h2h3h4_d <- as.data.frame(i_h2h3h4, row.names = NULL, col.names = NULL)
write.table(i_h2h3h4_d, "core_hydro_proteins.tsv", sep = "\t", row.names = FALSE)

```


```R
#core soil
#filtering soil samples 
soil<-filter_taxa(subset_samples(tom_all, substratum=="soil"), function (x) {sum(x > 0) > 0}, prune=TRUE)
soil
soil<-filter_taxa(subset_samples(soil, type=="sample"), function (x) {sum(x > 0) > 0}, prune=TRUE)
soil

#creating upset
entradaUpset <- otu_table(soil)
entradaUpset[entradaUpset>0] = 1
tomateUpset <- upset(as.data.frame(entradaUpset), nsets=40, nintersects=7, order.by = "degree", number.angles = 0, point.size = 2.5, line.size = 0.5, mainbar.y.label = "Genus intersections", 
                     sets.x.label = "Protein per sample",  text.scale = c(4,4,4,4,4,4))
tomateUpset
```


```R
#Get the core soil list
#Subset samples
soil1<-filter_taxa(subset_samples(soil, sample=="soil_1"), function (x) {sum(x > 0) > 0}, prune=TRUE)
soil2<-filter_taxa(subset_samples(soil, sample=="soil_2"), function (x) {sum(x > 0) > 0}, prune=TRUE)
soil3<-filter_taxa(subset_samples(soil, sample=="soil_3"), function (x) {sum(x > 0) > 0}, prune=TRUE)
soil4<-filter_taxa(subset_samples(soil, sample=="soil_4"), function (x) {sum(x > 0) > 0}, prune=TRUE)
soil5<-filter_taxa(subset_samples(soil, sample=="soil_5"), function (x) {sum(x > 0) > 0}, prune=TRUE)

# Save proteins list
soil1_prot <-as.data.frame(tax_table(soil1))
soil1_prot_v<-as.vector(row.names(soil1_prot))
soil2_prot <-as.data.frame(tax_table(soil2))
soil2_prot_v<-as.vector(row.names(soil2_prot))
soil3_prot <-as.data.frame(tax_table(soil3))
soil3_prot_v<-as.vector(row.names(soil3_prot))
soil4_prot <-as.data.frame(tax_table(soil4))
soil4_prot_v<-as.vector(row.names(soil4_prot))
soil5_prot <-as.data.frame(tax_table(soil5))
soil5_prot_v<-as.vector(row.names(soil5_prot))

# Obtain the intersection of all proteins
h2h3h4_l<-list(soil1_prot_v,soil2_prot_v,soil3_prot_v,soil4_prot_v,soil5_prot_v)
i_h2h3h4<-Intersect(h2h3h4_l)
i_h2h3h4_d<-as.data.frame(i_h2h3h4,row.names = NULL, col.names = NULL) 
write.table(i_h2h3h4_d, "core_SOIL_proteins.tsv", sep = "\t", row.names = FALSE)

```


```R
#Venn diagram function
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

#Load tables
core_hydro <-read.table("core_hydro_proteins.tsv", header = T)
core_hydro_v<-as.vector(core_hydro$i_h2h3h4)
core_soil <-read.table("core_SOIL_proteins.tsv", header = T)
core_soil_v<-as.vector(core_soil$i_h2h3h4)
```


```R
#create Venn diagram for core hydro and core soil
x2 <- list(core_hydro_v, core_soil_v)

# Display Venn diagram with customized settings
venn_cores <- display_venn(
  x2, 
  category.names = c("Core hydro", "Core soil"),
  lwd = 3,         
  cex = 2,         
  lty = 'blank',   
  fill = c("steelblue","rosybrown4"),  
  cat.cex = 1.7 
)
venn_cores
```

# Plot abundance of nitrogen fixation-related proteins


```R
# load phyloseq object with refseq tax_glom
# Loading data in phyloseq
otu <- as.matrix(read.table("taxglom_protein_otu_table.tsv", header=T, row.names=1)) 
OTU = otu_table(otu, taxa_are_rows=T)
taximat = as.matrix(read.table("taxglom_protein_tax_table.tsv", header=T, row.names=1)) 
taxi=tax_table(taximat)
tom_all= phyloseq(OTU, taxi)
sample_names(tom_all)
data= read.table("metada.csv", header=T, row.names=1, sep=",")
sampledata = sample_data(data.frame(sample=data$sample,subid=data$subsam, type=data$type,
                                    ID=data$ID, substratum=data$substratum, 
                                    row.names=sample_names(tom_all)))
tom_rs_taxglom= phyloseq(OTU, sampledata, taxi)
tom_rs_taxglom

# Estimate relative abundance
tom_rs_taxglom_r <- transform_sample_counts(tom_rs_taxglom, function(x) x / sum(x))
tom_rs_taxglom_r

                   
```


```R
# load phyloseq object with seed annotation data
#Loading data in phyloseq
otu <- as.matrix(read.table("global_ss_seed_counts.tsv", header=T, row.names=1)) 
OTU = otu_table(otu, taxa_are_rows=T)
taximat = as.matrix(read.table("global_ss_seed_onto.tsv", header=T, row.names=1)) 
taxi=tax_table(taximat)
tom_all= phyloseq(OTU, taxi)
sample_names(tom_all)
data= read.table("metada_fen.csv", header=T, row.names=1, sep=",")
sampledata = sample_data(data.frame(sample=data$sample,subid=data$subsam, type=data$type,
                                    ID=data$ID, substratum=data$substratum, 
                                    row.names=sample_names(tom_all)))
tom_all_seed= phyloseq (OTU, sampledata, taxi)
tom_all_seed
```


```R
# subset nitrogen_metabolism
nit_met_alls <- subset_taxa(tom_all_seed, level1=="nitrogen_metabolism")
nit_met_alls
# subset nitrogen fixation
fix_ni <- subset_taxa(nit_met_alls, level3=="nitrogen_fixation")
# save taxonomy table and OTU table for the subset
tom_tax_fix_ni <- tax_table(fix_ni)
# save the table
write.table(tom_tax_fix_ni, "tom_tax_fix_ni.tsv", sep = "\t")
# subset in refseq originating from seed
p1_ids <- read.table("tom_tax_fix_ni.tsv", header=TRUE, row.names=NULL, 
                     stringsAsFactors = FALSE, check.names=FALSE)
p1 <- p1_ids$row.names
# subset p1
sub_p1nit_met_alls <- subset_taxa(tom_all_rsq, rownames(tax_table(tom_all_rsq)) %in% p1)
sub_p1nit_met_alls
tom_tax_fix_ni_refseq <- tax_table(sub_p1nit_met_alls)
write.table(tom_tax_fix_ni_refseq, "tom_tax_fix_ni_refseq.tsv", sep = "\t")

# extract protein names
tt_sub_p1 <- as.data.frame(tax_table(sub_p1nit_met_alls))
# create a list of proteins
p_sub <- tt_sub_p1$protein1

# subset on my_taxglom to obtain relative abundances
sub_p2 <- subset_taxa(tom_rs_taxglom_r, protein1 %in% p_sub)
sub_p2

# extract data for plotting
p <- plot_bar(sub_p2, "protein1")
yx <- p$data
yx <- as.data.frame(yx)
# sort by abundance
yx <- yx[order(-yx$Abundance),]

# replace NA with zeros
yx[is.na(yx)] <- 0
# remove hypothetical proteins
yx2 <- yx[-grep("hypothetical_protein", yx$protein1),]
write.table(yx2, file='yx2.tsv', quote=FALSE, sep='\t', col.names = NA)

# load table
yx2 <- read.table("yx2.tsv", header=TRUE, row.names=NULL, sep = "\t")
# clean table
yx3 <- yx2[yx2$Abundance != 0,]

# order samples
orden_samples <- factor(yx3$Sample,levels= c("soil5ss1", "soil5ss2", "soil5ss3","soil5",
                                             "soil4ss1", "soil4ss2", "soil4ss3","soil4", 
                                             "soil3ss1", "soil3ss2", "soil3ss3","soil3", 
                                             "soil2ss1", "soil2ss2", "soil2ss3", "soil2",
                                             "soil1ss1", "soil1ss2", "soil1ss3","soil1", 
                                             "hydro4ss1", "hydro4ss2", "hydro4ss3","hydro4",
                                             "hydro3ss1", "hydro3ss2", "hydro3ss3", "hydro3",
                                             "hydro2ss1", "hydro2ss2", "hydro2ss3","hydro2"))

# Plot relative abundance of phyla
ggplot(yx3, aes(x = reorder(protein1, desc(Abundance)), 
                y = orden_samples, 
                fill = Abundance)) + 
geom_raster() + 
theme_minimal() + 
theme(text = element_text(size = 14)) +
scale_fill_gradient(low = "#e4e3e5ff", high = "#604a7bff", trans = "log10") +
theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

# Plot the taxonomic assignment of proteins related to nitrogen fixation


```R
# load tables into R 
table_1 <- read.table("nit_fix_seq_name_clean_v2.txt", header=FALSE, row.names=NULL, check.names=FALSE)
colnames(table_1) <- c("id", "protein", "contig_a")
table_2 <- read.table("nit_fix_ann.tsv", header=FALSE, row.names=NULL, check.names=FALSE)
colnames(table_2) <- c("contig_b", "tax", "tax_id")
#nit_fix_seq_name_clean_v2.txt contains M5nr annotation of proteins related to nitrogen fixation
#nit_fix_ann.tsv contains taxonomic annotation of proteins related to nitrogen fixation

# merge the tables
m1 <- cbind(table_1, table_2)

# load taxonomic table from refseq for nitrogen fixation proteins
nitfix_tax_rs <- read.table("tom_tax_fix_ni_refseq.tsv", header=TRUE, row.names=NULL, sep="\t", check.names=FALSE)

# merge taxonomic data with annotation data
m2 <- merge(nitfix_tax_rs, m1, by="id")

# calculate the occurrence and relative occurrence of each protein in different taxa
ctp <- m2 %>% group_by(protein1, tax) %>% summarize(count = n()) %>% mutate(relative_occurrence = count/sum(count))

# heatmap
# remove hypothetical proteins
ctp2 <- ctp[-grep("hypothetical_protein", ctp$protein1),]

# order proteins
orden <- factor(ctp2$protein1, levels = c("putative_nitrogen_fixation_protein_NifT", "nitrogen_fixation_protein_NifX", "ferredoxin_III_nif_specific", "nitrogenase_stabilizing_protective_protein_NifW", "nitrogen_fixation_protein_NifZ", "nitrogen_fixation_protein_NifQ", "NifX_associated_nitrogen_fixation_protein", "nitrogenase_iron_protein", "nitrogenase_molybdenum_iron_protein_alpha_chain", "homocitrate_synthase", "nitrogenase_cofactor_biosynthesis_protein_NifB", "nitrogenase_iron_molybdenum_cofactor_biosynthesis_protein_NifN", "nitrogenase_iron_molybdenum_cofactor_biosynthesis_protein_NifE", "nitrogenase_molybdenum_iron_protein_subunit_beta", "nif_specific_transcriptional_activator_NifA", "iron_sulfur_cluster_assembly_accessory_protein", "4Fe_4S_binding_protein"))

# order taxa by total relative occurrence
ctp2_summed <- ctp2 %>%
  group_by(tax) %>%
  summarise(total_relative_occurrence = sum(relative_occurrence)) %>%
  arrange(desc(total_relative_occurrence))

# reorder taxa for plotting
ctp2_ordered <- ctp2 %>%
  mutate(tax = factor(tax, levels = rev(ctp2_summed$tax)))

# plot the heatmap
ggplot(ctp2_ordered, aes(x = factor(orden, levels = rev(levels(orden))), 
                          y = tax,
                          fill = relative_occurrence)) + 
  geom_raster() + 
  theme_minimal() + 
  theme(text = element_text(size = 14)) +
  scale_fill_gradient(low = "#e4ede6ff", high = "#5a8862ff", na.value = "black", trans = "log10") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
```
