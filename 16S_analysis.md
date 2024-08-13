```R
#files 
#allhydro_FINAL_c1m1b1.otu.tsv
#taxonomy_FINAL_c1m1b1.tsv
#metadata_hydro_order.tsv
```


```R
#load packages
library(ggplot2)
library(vegan)
library(phyloseq)
library(RColorBrewer)
library(reshape2)
library(ggsignif)
library(ggplot2)
library(ape)
library(VennDiagram)
library(tidyverse)
library(tibble)
library(ggnewscale)
library(forcats)
library(scales)
library(UpSetR)
library(grDevices)
library(gridExtra)
```

# 16S rRNA Amplicon Analysis


```R
#Load phyloseq object
otu <- as.matrix(read.table("allhydro_FINAL_c1m1b1.otu.tsv", header=T, row.names=1)) 
OTU = otu_table(otu, taxa_are_rows=T)
taximat = as.matrix(read.table("taxonomy_FINAL_c1m1b1.tsv", header=T, row.names=1, sep="\t")) 
taxi=tax_table(taximat)
tom_hydro= phyloseq(OTU, taxi)
sample_names(tom_hydro)
data= read.table("metadata_hydro_order.tsv", header=T, row.names=1, sep="\t")
head(data)
sampledata = sample_data(data.frame(sample=data$sample,sam_name=data$sam_name, Rh_En=data$R_E,
                                    treatment=data$treatment,tag=data$tag,
                                    treatment_c=data$treatment_c, type=data$type,
                                    pH=data$pH, stem_length= data$stem_length,
                                    wet_weight=data$wet_weight, dry_weight=data$dry_weight,
                                    foliar_surface=data$foliar_surface, stem_diameter=data$stem_diameter,
                                    row.names=sample_names(tom_hydro)))

tom_hydro= phyloseq (OTU, sampledata, taxi)
print("objeto de phyloseq cargado")
tom_hydro
```

# Estimation of alpha diversity


```R
#Calculation of diversity indices
div_table<-estimate_richness(tom_hydro , measures=c("Observed", "Chao1", "Shannon", "Simpson"))
write.table(div_table, "index_diversidad_hydro.tsv", sep="\t")
#Obtain data from plot
p <- plot_richness(tom_hydro, x="treatment_c", color="treatment_c", 
                        measures=c("Observed", "Chao1", "Shannon", "Simpson"))
yx <- p$data
yx <- as.data.frame(yx)

#Extract values for each diversity index
shannon <- yx %>% filter(variable == "Shannon")
```


```R
#Shannon plot
#Prepare the data and create the boxplot
orden <- factor(shannon$treatment_c, levels = c('f50', 'f100', 'nc', 'treatment'))
cajas_dob <- ggplot(data = shannon, aes(x = orden, y = value, colour = factor(treatment_c))) +
  geom_boxplot() +
  geom_point(size = 3) +
  scale_color_manual(values = c('gray53', 'coral3', 'palegreen4', 'darkgoldenrod2')) +theme_light() 

#Add significance lines using geom_signif
cajas_dob + geom_signif(test = 't.test',comparisons = list(c("f50", "f100")), 
                                           map_signif_level = TRUE, textsize = 6, y_position = 6) +
                         geom_signif(test = 't.test', comparisons = list(c("f50", "treatment")), 
                                           map_signif_level = TRUE, textsize = 6, y_position = 6.75) +
                                      geom_signif(test = 't.test',comparisons = list(c("f100", "treatment")), 
                                           map_signif_level = TRUE, textsize = 6, y_position = 6.25) +
                                                  geom_signif(test = 't.test',comparisons = list(c("nc", "treatment")), 
                                           map_signif_level = TRUE, textsize = 6, y_position = 6)+
theme(text = element_text(size=30)) + 
theme(legend.position = "none") + labs(title = "Shannon", y = "Shannon", x = NULL)
```

# Plotting Dry Weight


```R
# Load dry weight data
peso_seco= read.table("peso_seco.csv", header=T, row.names=NULL, sep=",")
head(peso_seco)
```


```R
#Prepare the data and generate the boxplot
orden <- factor(peso_seco$treatment_c, levels = c('f50', 'f100', 'nc', 'treatment'))
cajas_dob <- ggplot(data = peso_seco, aes(x = orden, y = dry_weight, colour = factor(treatment_c))) +
  geom_boxplot() +
geom_point(size = 5) +
  scale_color_manual(values = c('gray53', 'coral3', 'palegreen4', 'darkgoldenrod2')) +theme_light() 

#Add significance lines using geom_signif
cajas_dob + geom_signif(test = 't.test',comparisons = list(c("f50", "f100")), 
                                           map_signif_level = TRUE, textsize = 6, y_position = 6) +
                         geom_signif(test = 't.test',comparisons = list(c("f50", "treatment")), 
                                           map_signif_level = TRUE, textsize = 6, y_position = 7.2) +
                                      geom_signif(test = 't.test',comparisons = list(c("f100", "treatment")), 
                                           map_signif_level = TRUE, textsize = 6, y_position = 6.6) +
                                                  geom_signif(test = 't.test',comparisons = list(c("nc", "treatment")), 
                                           map_signif_level = TRUE, textsize = 6, y_position = 6) +
theme(text = element_text(size=30)) + 
theme(legend.position = "none") + labs(title = "Dry weight", y = "Dry weight", x = NULL)
```

# CAP Ordination Plot


```R
# Perform Canonical Analysis of Principal coordinates (CAP) ordination
# using Bray-Curtis dissimilarity matrix and specified environmental variables
fun_CAP_bray = ordinate(tom_hydro, "CAP", "bray", ~dry_weight +
                      stem_length + foliar_surface + stem_diameter + pH)

# Create a plot of the CAP ordination, using the specified shape and color
# aesthetics for the samples
cap_plot = plot_ordination(tom_hydro, shape = "Rh_En", fun_CAP_bray, type = "samples", color = "treatment_c") + 
  scale_color_manual(values = c('gray53', 'coral3', 'palegreen4', 'darkgoldenrod2')) + theme_bw() +
  geom_point(size = 7) + theme(text = element_text(size = 20))

# Extract the scores for the biplot axes
arrowmat <- vegan::scores(fun_CAP_bray, display = "bp")
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the aesthetic mapping for arrows
arrow_map <- aes(xend = CAP1,
                 yend = CAP2,
                 x = 0,
                 y = 0,
                 shape = NULL,
                 color = black,
                 label = labels)

# Define the aesthetic mapping for labels
label_map <- aes(x = 1.3 * CAP1,
                 y = 1.3 * CAP2,
                 shape = NULL,
                 color = NULL,
                 label = labels)
# Add arrows and labels to the CAP plot
arrowhead = arrow(length = unit(0.02, "npc"))

capl_plot <- cap_plot  +
  geom_segment(
    mapping = arrow_map,
    size = .8,
    data = arrowdf,
    color = "black",
    arrow = arrowhead
  ) +
  geom_text(
    mapping = label_map,
    size = 6, 
    data = arrowdf,
    show.legend = FALSE
  )

# Display the final CAP plot with arrows and labels
capl_plot

# Perform permutation-based ANOVA to test the significance of the CAP axes
anx = anova(fun_CAP_bray, permutations = 9999, alpha = 0.05)
anx

# Perform permutation-based ANOVA for each term individually
any = anova(fun_CAP_bray, by = 'terms', permutations = 9999, alpha = 0.05)
any
```

# ADONIS calculation


```R
# Convert the OTU table to a matrix
tom_ott <- as.matrix(otu_table(tom_hydro))
tom_ott_t <- t(tom_ott)
dist_matrix <- vegdist(tom_ott_t, method = "bray")
# Format the distance matrix as a matrix
dist_matrix_matrix <- as.matrix(dist_matrix)

# Filter the sample data for treatment and nc groups
sample_data_tom <- data.frame(sample_data(tom_hydro))
treatment_nc <- sample_data_tom %>% filter(treatment_c %in% c("treatment", "nc"))
# Get the sample names for the group
treatment_nc_samples <- rownames(treatment_nc)
# Filter the distance matrix for the samples names
filtered_dist_matrix <- dist_matrix_matrix[treatment_nc_samples, treatment_nc_samples]
# Filter the sample data for the common samples
treatment_data_filtered <- treatment_nc[treatment_nc_samples, ]
# Perform ADONIS analysis within the "treatment" group
adonis_result_nc <- adonis2(filtered_dist_matrix ~ treatment, data = treatment_data_filtered, permutations = 9999)
adonis_result_nc

# Filter the sample data for treatment and f50 groups
sample_data_tom <- data.frame(sample_data(tom_hydro))
treatment_f50 <- sample_data_tom %>% filter(treatment_c %in% c("treatment", "f50"))
# Get the sample names for the group
treatment_f50_samples <- rownames(treatment_f50)
# Filter the distaf50e matrix for the samples names
filtered_dist_matrix <- dist_matrix_matrix[treatment_f50_samples, treatment_f50_samples]
# Filter the sample data for the common samples
treatment_data_filtered <- treatment_f50[treatment_f50_samples, ]
# Perform ADONIS analysis within the "treatment" group
adonis_result_f50 <- adonis2(filtered_dist_matrix ~ treatment, data = treatment_data_filtered, permutations = 9999)
adonis_result_f50

# Filter the sample data for treatment and f100 groups
sample_data_tom <- data.frame(sample_data(tom_hydro))
treatment_f100 <- sample_data_tom %>% filter(treatment_c %in% c("treatment", "f100"))
# Get the sample names for the group
treatment_f100_samples <- rownames(treatment_f100)
# Filter the distaf100e matrix for the samples names
filtered_dist_matrix <- dist_matrix_matrix[treatment_f100_samples, treatment_f100_samples]
# Filter the sample data for the common samples
treatment_data_filtered <- treatment_f100[treatment_f100_samples, ]
# Perform ADONIS analysis within the "treatment" group
adonis_result_f100 <- adonis2(filtered_dist_matrix ~ treatment, data = treatment_data_filtered, permutations = 9999)
adonis_result_f100

```

# Taxomony composition


```R
#Phylum
```


```R
#Most abundant phylum
#Tax glom at Phylum level
phy_tom <- tax_glom(tom_hydro, "Phylum")
phy_tom

#Estimate relative abundance
rel_phy <- transform_sample_counts(phy_tom, function(x) x / sum(x))
rel_phy
```


```R
# Get the file with taxa names and taxa counts
# Extract the OTU table and taxonomy table from the phyloseq object
tom_otu_phy <- otu_table(rel_phy)
tom_tax_phy <- tax_table(rel_phy)

# Combine OTU table and taxonomy table
tom_otu_tax <- cbind(tom_otu_phy, tom_tax_phy)
# Select columns from 1 to 34
tom_otu_tax <- tom_otu_tax[, 1:34]

# Write the combined OTU and taxonomy table to a TSV file
write.table(tom_otu_tax, "tom_otu_tax_phylum.tsv", sep = "\t")

# Read the TSV file back into R
tom_otu_tax <- read.table("tom_otu_tax_phylum.tsv", header = TRUE, row.names = 1, 
                                 stringsAsFactors = FALSE, check.names = FALSE)

# Create a copy of the data frame
tom_otu_tax_b <- tom_otu_tax

# Melt the data frame to long format for specific samples
m_tom_otu_tax <- melt(tom_otu_tax, measure.vars = c(  
     "hydro1", "hydro2", "hydro3", "hydro4", "hydro5",  "hydro6", "hydro7", "hydro8", "hydro13", "hydro14", 
    "hydro15", "hydro16", "hydro21", "hydro22", "hydro23", "hydro24", "hydro29", "hydro30", "hydro31",
    "hydro32", "hydro9",  "hydro10", "hydro11", "hydro12", "hydro17", "hydro18", "hydro19", "hydro20", "hydro25",
    "hydro26", "hydro27", "hydro28"))

# Rename columns for clarity
colnames(m_tom_otu_tax) <- c("Kingdom", "Phylum", "Sample", "Relative_abundance")

# Select columns 2 to 4
m_tom_otu_tax <- m_tom_otu_tax[, 2:4]

# Write the melted data frame to a TSV file
write.table(m_tom_otu_tax, "tax_phylum_melt.tsv", sep = "\t")
```


```R
# Calculate the sum of Relative_abundance for each Phylum
datos_colapsados <- m_tom_otu_tax %>%
  group_by(Phylum) %>%
  summarise(sum_abundance = sum(Relative_abundance))

# Identify Phyla with a sum of abundance less than 0.001
phyla_low <- datos_colapsados %>%
  filter(sum_abundance < 0.001) %>%
  pull(Phylum)

# Create a new column "Category" in m_tom_otu_tax to classify low-abundance Phyla
m_tom_otu_tax <- m_tom_otu_tax %>%
  mutate(Category = ifelse(Phylum %in% phyla_low, "Low_abundance", Phylum))

# Collapse the data by Sample and Category, summing the Relative_abundance
collapsed_data <- m_tom_otu_tax %>%
  group_by(Sample, Category) %>%
  summarise(sum_abundance = sum(Relative_abundance))

# Define the order of samples for the heatmap
order <- factor(collapsed_data$Sample, levels = c("hydro5", "hydro6", "hydro7", "hydro8", 
                                    "hydro1", "hydro2", "hydro3", "hydro4",
                                    "hydro13", "hydro14", "hydro15", "hydro16", 
                                    "hydro21", "hydro22", "hydro23", "hydro24",
                                    "hydro29", "hydro30", "hydro31", "hydro32",
                                    "hydro9", "hydro10", "hydro11", "hydro12", "hydro17", 
                                    "hydro18", "hydro19", "hydro20", "hydro25", 
                                    "hydro26", "hydro27", "hydro28"))

# Create a heatmap
heatmap <- ggplot(collapsed_data, aes(x = order, y = reorder(Category, sum_abundance), fill = sum_abundance)) +
  geom_raster() +
  theme_minimal() +
  theme(text = element_text(size = 30), axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient(low = "#bebebe", high = "#000000", na.value = "white", trans = "log10")
heatmap
```

# Correlogram of genus abundance with plant phenotype


```R
# Agglomerate OTUs based on genus
tom_hydro_genus_ms <- tax_glom(tom_hydro, taxrank = "Genus")
tom_hydro_genus_ms

# Estimate relative abundance
tom_hydro_genus_rel <- transform_sample_counts(tom_hydro_genus_ms, function(x) x / sum(x))
tom_hydro_genus_rel
```


```R
# Save taxonomy table
tom_hydro_tax <- tax_table(tom_hydro_genus_ms)
write.table(tom_hydro_tax, "tom_hydro_tax.tsv", sep = "\t")
# Calculate relative abundances
tom_hydro_genus_rel <- transform_sample_counts(tom_hydro_genus_ms, function(x) x / sum(x))
tom_hydro_genus_rel
tom_hydro_otu_table <- otu_table(tom_hydro_genus_rel)
write.table(tom_hydro_otu_table, "tom_hydro_otu_table.tsv", sep = "\t")
# Load tables
hydro_tax <- as.matrix(read.table("tom_hydro_tax.tsv", header = TRUE, row.names = 1))
hydro_otu_table <- as.matrix(read.table("tom_hydro_otu_table.tsv", header = TRUE, row.names = 1))
# Merge taxonomy and OTU count tables based on the row names containing the OTU numbers
m <- merge(hydro_tax, hydro_otu_table, by = 'row.names', all = TRUE)
# Combine OTU taxonomy into a single column
m$tax_g <- paste(m$Phylum, "|", m$Class, "|", m$Order, "|", m$Family, "|", m$Genus)
# Transpose the table
mt <- t(m)
# Remove taxonomy rows to retain combined taxonomy
mtR <- mt[-c(2:6), ]
# Remove OTU number row
mtRnotu <- mtR[-c(1, 2), ]
# Save combined taxonomy names
genera_names <- as.character(mtRnotu[33, ])
# Add taxonomy as column headers
colnames(mtRnotu) <- genera_names
# Remove genus row as it's now in the header
mtR2 <- mtRnotu[-c(33), ]
# Save relative abundance table
write.table(mtR2, "relative_abundances_for_corr.tsv", sep = "\t")
# Load phenotype data
fenotipo0 <- sample_data(tom_hydro)
# Remove non-numeric columns
fenotipo <- fenotipo0[, c(7:12)]
# Merge phenotype data with relative abundance table by sample number
fen_rel <- merge(mtR2, fenotipo, by = 'row.names', all = TRUE)
# Remove row names column
fen_rel$Row.names <- NULL
# Save the combined phenotype and abundance data
write.csv(fen_rel, "fen_rel.csv")                                           
```


```R
#Correlations were calculated using Spearman's method with Benjamini-Hochberg correction for false discovery rate

col <- colorRampPalette(c("deepskyblue4", "white", "seagreen3"))(20)
rmatrix <- read.csv('fen_rel.csv', header=T, row.names=1, sep=",")
rcorr <- corr.test(rmatrix, method = "spearman", adjust = "BH", alpha = 0.05)
dfr2<-rcorr$r

write.csv(dfr2, "rcorr2.csv", row.names = TRUE)
```


```R
# Load the correlation table
corr <- read.table("rcorr2.csv", header = TRUE, row.names = 1, sep = ",")
# Remove columns containing genera
corr_clean <- corr[, c(819:824)]
# Remove rows containing phenotype data
corr_clean2 <- corr_clean[1:818, ]
# Convert row names to a column and clean the table
corr_clean2_nrn <- tibble::rownames_to_column(corr_clean2, "genus")
corr_clean_c <- corr_clean2_nrn[, c(1, 3:7)]
# Convert the data to a format suitable for the heatmap
corr_clean_melt <- melt(corr_clean_c, id.vars = c("genus"))
# Load and reverse the order of correlation values
order_corrval <- read.table("corr_order2.txt", sep = "\t")
order_corrval <- as.character(order_corrval$V1)
order_corrval_r <- rev(order_corrval)
order_y <- factor(corr_clean_melt$genus, levels = order_corrval_r)
# Plot the correlation values
my_palette <- c("#E1BE6A", "white", "#40B0A6")

ggplot(corr_clean_melt, aes(x = variable, y = order_y, fill = value)) + 
  geom_raster() + 
  theme_minimal() + 
  theme(text = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradientn(colours = my_palette, values = scales::rescale(c(-0.7, 0, 0.7)))

```

# Plot the abundance of all genera, grouped by control and treatment


```R
# Merge samples by treatment category
tom_hydro_merge <- merge_samples(tom_hydro, "treatment_c")
tom_hydro_merge

# Agglomerate OTUs based on genus
tom_hydro_merge_genus <- tax_glom(tom_hydro_merge, taxrank = "Genus")
tom_hydro_merge_genus

# Estimate relative abundance
tom_hydro_merge_genus_r <- transform_sample_counts(tom_hydro_merge_genus, function(x) x / sum(x))
tom_hydro_merge_genus_r

# Extract data for plotting
p <- plot_bar(tom_hydro_merge_genus_r, "Genus")
yx <- p$data
yx <- as.data.frame(yx)
yx <- yx[order(-yx$Abundance),]

# Combine OTU taxonomy into a single column
yx$tax_g <- paste(yx$Phylum, "|", yx$Class, "|", yx$Order, "|", yx$Family, "|", yx$Genus)

# Define the order of samples
order <- factor(yx$Sample, levels = c("f50", "f100", "nc", "treatment"))

# Load and reverse the order of correlation values
order_corrval <- read.table("corr_order2.txt", sep = "\t")
order_corrval <- as.character(order_corrval$V1)
order_corrval_r <- rev(order_corrval)

# Define the order of genera for the y-axis
order_y <- factor(yx$tax_g, levels = order_corrval_r)

# Plot relative abundance of genera
ggplot(yx, aes(x = order, y = order_y, fill = Abundance)) + 
  geom_raster() + 
  theme_minimal() + 
  theme(text = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient(low = "#bebebe", high = "#000000", na.value = "white", trans = "log10")

```

# Plot of genera with low abundance collapsed into a single category


```R
# Get the file with taxa names and taxa counts
tom_otu_gen <- otu_table(tom_hydro_genus_rel)
tom_tax_gen <- tax_table(tom_hydro_genus_rel)
tom_otu_tax <- cbind(tom_otu_gen, tom_tax_gen)
tom_otu_tax <- tom_otu_tax[, c(1:32, 38)]
write.table(tom_otu_tax, "tom_otu_tax_genus.tsv", sep = "\t")

# Read the file and reshape the data for analysis
tom_otu_tax <- read.table("tom_otu_tax_genus.tsv", header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
m_tom_otu_tax <- melt(tom_otu_tax, measure.vars = c(
  "hydro1", "hydro2", "hydro3", "hydro4", "hydro5", "hydro6", "hydro7", "hydro8", "hydro13", "hydro14", 
  "hydro15", "hydro16", "hydro21", "hydro22", "hydro23", "hydro24", "hydro29", "hydro30", "hydro31",
  "hydro32", "hydro9", "hydro10", "hydro11", "hydro12", "hydro17", "hydro18", "hydro19", "hydro20", 
  "hydro25", "hydro26", "hydro27", "hydro28"
))
colnames(m_tom_otu_tax) <- c("Genus", "Sample", "Relative_abundance")

# Calculate the sum of Relative_abundance by Genus
datos_colapsados <- m_tom_otu_tax %>%
  group_by(Genus) %>%
  summarise(sum_abundance = sum(Relative_abundance))

# Identify the genera with a sum less than 0.05
genus_low <- datos_colapsados %>%
  filter(sum_abundance < 0.05) %>%
  pull(Genus)

# Create a new column "Category" in m_tom_otu_tax
low_m_tom_otu_tax <- m_tom_otu_tax %>%
  mutate(Category = ifelse(Genus %in% genus_low, "Low abundance", Genus))
low_m_tom_otu_tax <- low_m_tom_otu_tax[!(low_m_tom_otu_tax$Genus == "_g__"), ]
                                               
# Collapse the data by Sample and Category, summing Relative_abundance
collapsed_data <- low_m_tom_otu_tax %>%
  group_by(Sample, Category) %>%
  summarise(sum_abundance = sum(Relative_abundance))

# Order the samples
order <- factor(collapsed_data$Sample, levels = c(
  "hydro5", "hydro6", "hydro7", "hydro8", "hydro1", "hydro2", "hydro3", "hydro4",
  "hydro13", "hydro14", "hydro15", "hydro16", "hydro21", "hydro22", "hydro23", "hydro24",
  "hydro29", "hydro30", "hydro31", "hydro32", "hydro9", "hydro10", "hydro11", "hydro12",
  "hydro17", "hydro18", "hydro19", "hydro20", "hydro25", "hydro26", "hydro27", "hydro28"
))

# Create heatmap
heatmap <- ggplot(collapsed_data, aes(x = order, y = reorder(Category, sum_abundance), fill = sum_abundance)) +
  geom_raster() +
  theme_minimal() +
  theme(text = element_text(size = 30), axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient(low = "#bebebe", high = "#000000", na.value = "white", trans = "log10")
heatmap
```

# Venn diagram with genera


```R
#funcion para mostrar el diagrama
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
```


```R
# Get the lists of genera in each type of control and treatment

# Subset of fertilization controls (all controls)
# This subset is for hydroponics, considering all fertilization controls
genus_cf <- filter_taxa(subset_samples(tom_hydro_genus_ms, tag == "control_f"), function(x) {sum(x > 0) > 0}, prune = TRUE)
genus_cf

# Subset of fertilization controls (50%)
# This subset is for hydroponics, considering fertilization controls with 50% fertilizer
genus_f50 <- filter_taxa(subset_samples(tom_hydro_genus_ms, treatment_c == "f50"), function(x) {sum(x > 0) > 0}, prune = TRUE)
genus_f50

# Subset of fertilization controls (100%)
# This subset is for hydroponics, considering fertilization controls with 100% fertilizer
genus_f100 <- filter_taxa(subset_samples(tom_hydro_genus_ms, treatment_c == "f100"), function(x) {sum(x > 0) > 0}, prune = TRUE)
genus_f100

# Subset of sterility controls
# This subset is for hydroponics, considering sterility controls
genus_nc <- filter_taxa(subset_samples(tom_hydro_genus_ms, treatment_c == "nc"), function(x) {sum(x > 0) > 0}, prune = TRUE)
genus_nc

# Subset of hydroponic treatments
# This subset is for hydroponics, considering the inoculants
genus_treatment <- filter_taxa(subset_samples(tom_hydro_genus_ms, treatment_c == "treatment"), function(x) {sum(x > 0) > 0}, prune = TRUE)
genus_treatment

# Save the genera in vectors
genus_treatment_tax <- as.data.frame(tax_table(genus_treatment))
genus_treatment_lab <- apply(genus_treatment_tax, 1, function(x) paste(x[-1], collapse = " | "))
genus_treatment_v <- as.vector(genus_treatment_lab)
write.table(genus_treatment_lab, "genera_treatment.tsv", sep = "\t", row.names = FALSE)

# Save the genera in vectors
genus_f100_tax <- as.data.frame(tax_table(genus_f100))
genus_f100_lab <- apply(genus_f100_tax, 1, function(x) paste(x[-1], collapse = " | "))
genus_f100_lab_v <- as.vector(genus_f100_lab)
write.table(genus_f100_lab, "genera_f100.tsv", sep = "\t", row.names = FALSE)

# Save the genera in vectors
genus_f50_tax <- as.data.frame(tax_table(genus_f50))
genus_f50_lab <- apply(genus_f50_tax, 1, function(x) paste(x[-1], collapse = " | "))
genus_f50_lab_v <- as.vector(genus_f50_lab)
write.table(genus_f50_lab, "genera_f50.tsv", sep = "\t", row.names = FALSE)

# Save the genera in vectors
genus_nc_tax <- as.data.frame(tax_table(genus_nc))
genus_nc_lab <- apply(genus_nc_tax, 1, function(x) paste(x[-1], collapse = " | "))
genus_nc_lab_v <- as.vector(genus_nc_lab)
write.table(genus_nc_lab, "genera_nc.tsv", sep = "\t", row.names = FALSE)

# Create Venn diagram
x2 <- list(genus_treatment_v, genus_nc_lab_v, genus_f50_lab_v, genus_f100_lab_v)

venn_hydro <- display_venn(x2, category.names = c("treatment", "NC", "f50", "f100"),
                           lwd = 3, cex = 3,
                           lty = 'blank',
                           fill = c("darkgoldenrod", "palegreen4", "coral3", "cornsilk2"))
```

# Calculation of overrepresented genera in treatments compared to controls


```R
# filtering controls
controls <- filter_taxa(subset_samples(tom_hydro_genus_ms, 
                                  type=="control"), function (x) {sum(x > 0) > 0}, prune=TRUE)
controls

# here, extract the row names of OTUs
# save the number of otus
treat_otus <- as.vector(row.names(tax_table(treat)))
controls_otus <- as.vector(row.names(tax_table(controls)))

# creating the lists for each set
l_t_allcon <- list(treat_otus, controls_otus)

# obtain the intersection of controls and treatments
tyconall <- Intersect(l_t_allcon)
tyconall
tyconall_d <- as.data.frame(tyconall, row.names = NULL)
nrow(tyconall_d)

# obtain the OTUs shared by controls and treatments
subset_comp_tyconall <- subset_taxa(tom_hydro_genus_ms, rownames(tax_table(tom_hydro_genus_ms)) %in% tyconall)
subset_comp_tyconall

# enrichment analysis of controls and treatments
alpha = 0.05

merge_phyloseq(subset_samples(subset_comp_tyconall, type=="treatment"), 
               subset_samples(subset_comp_tyconall, type=="control")) -> sar.shu.phy

sar.shu.ds <- phyloseq_to_deseq2(sar.shu.phy, ~tag)
sar.shu.ds <- DESeq(sar.shu.ds, test="Wald", fitType="local")
sar.shu.ds.res <- results(sar.shu.ds, cooksCutoff = FALSE)
# write.csv(sar.shu.ds.res, "sar.shu.ds.res_tratamientosvscesterilidad.csv" )
sigtab.sar.shu <- sar.shu.ds.res[which(sar.shu.ds.res$padj < alpha), ]
sigtab.sar.shu <- cbind(as(sigtab.sar.shu, "data.frame"), 
                        as(tax_table(sar.shu.phy)[rownames(sigtab.sar.shu), ], "matrix"))
sigtab.sar.shu.x <- tapply(sigtab.sar.shu$log2FoldChange, sigtab.sar.shu$Genus, 
                           function(x) max(x))
sigtab.sar.shu.x <- sort(sigtab.sar.shu.x, TRUE)
sigtab.sar.shu$Genus <- factor(as.character(sigtab.sar.shu$Genus), 
                               levels=names(sigtab.sar.shu.x))
write.csv(sigtab.sar.shu, "log2_tyallcontrols_05.csv")
```

# Core hydroponics calculation


```R
#filtering treatment 
treatment_genera<-filter_taxa(subset_samples(tom_hydro_genus_ms, treatment_c=="treatment"), function (x) {sum(x > 0) > 0}, prune=TRUE)
treatment_genera

#creating upset
entradaUpset <- otu_table(treatment_genera)
entradaUpset[entradaUpset>0] = 1
tomateUpset <- upset(as.data.frame(entradaUpset), nsets=40, nintersects=6, order.by = "degree", number.angles = 0, point.size = 2.5, line.size = 0.5, mainbar.y.label = "Genus intersections", 
                     sets.x.label = "Genus per sample",  text.scale = c(4,4,4,4,4,4))
tomateUpset
```


```R
#defining functions to obtain the list of genera in the core

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
```


```R
#get list of genus present in treatments
#Subset treatment hydro2_1
hydro2_1_e<-filter_taxa(subset_samples(tom_hydro_genus_ms, sample=="hydro_s2-2_e"), function (x) {sum(x > 0) > 0}, prune=TRUE)
#Subset treatment hydro2_2
hydro2_2_e<-filter_taxa(subset_samples(tom_hydro_genus_ms, sample=="hydro_s2-3_e"), function (x) {sum(x > 0) > 0}, prune=TRUE)
#Subset treatment hydro3_1
hydro3_1_e<-filter_taxa(subset_samples(tom_hydro_genus_ms, sample=="hydro_s3-2_e"), function (x) {sum(x > 0) > 0}, prune=TRUE)
#Subset treatment hydro3_2
hydro3_2_e<-filter_taxa(subset_samples(tom_hydro_genus_ms, sample=="hydro_s3-3_e"), function (x) {sum(x > 0) > 0}, prune=TRUE)
#Subset treatment hydro4_1
hydro4_1_e<-filter_taxa(subset_samples(tom_hydro_genus_ms, sample=="hydro_s4-1_e"), function (x) {sum(x > 0) > 0}, prune=TRUE)
#Subset treatment hydro4_2
hydro4_2_e<-filter_taxa(subset_samples(tom_hydro_genus_ms, sample=="hydro_s4-2_e"), function (x) {sum(x > 0) > 0}, prune=TRUE)
#Subset treatment hydro2_1
hydro2_1_r<-filter_taxa(subset_samples(tom_hydro_genus_ms, sample=="hydro_s2-2_r"), function (x) {sum(x > 0) > 0}, prune=TRUE)
#Subset treatment hydro2_2
hydro2_2_r<-filter_taxa(subset_samples(tom_hydro_genus_ms, sample=="hydro_s2-3_r"), function (x) {sum(x > 0) > 0}, prune=TRUE)
#Subset treatment hydro3_1
hydro3_1_r<-filter_taxa(subset_samples(tom_hydro_genus_ms, sample=="hydro_s3-2_r"), function (x) {sum(x > 0) > 0}, prune=TRUE)
#Subset treatment hydro3_2
hydro3_2_r<-filter_taxa(subset_samples(tom_hydro_genus_ms, sample=="hydro_s3-3_r"), function (x) {sum(x > 0) > 0}, prune=TRUE)
#Subset treatment hydro4_1
hydro4_1_r<-filter_taxa(subset_samples(tom_hydro_genus_ms, sample=="hydro_s4-1_r"), function (x) {sum(x > 0) > 0}, prune=TRUE)
#Subset treatment hydro4_2
hydro4_2_r<-filter_taxa(subset_samples(tom_hydro_genus_ms, sample=="hydro_s4-2_r"), function (x) {sum(x > 0) > 0}, prune=TRUE)

```


```R
# List of sample names
samples <- c("hydro2_1_e", "hydro2_2_e", "hydro3_1_e", "hydro3_2_e", "hydro4_1_e", "hydro4_2_e", 
             "hydro2_1_r", "hydro2_2_r", "hydro3_1_r", "hydro3_2_r", "hydro4_1_r", "hydro4_2_r")

# Initialize an empty list to store the resulting vectors
treatments_genera_list <- list()

# Initialize an empty list to store each tax_lab
tax_lab_list <- list()

# Loop to process each sample
for (sample in samples) {
  # Convert the tax table to a data frame
  tax_df <- as.data.frame(tax_table(get(sample)))
  
  # Extract row names corresponding to otu numbers
  tax_vec <- rownames(tax_df)
  
  # Add the resulting vector to the list
  treatments_genera_list[[sample]] <- tax_vec
    
   # Create a dynamic variable name and store the vector
  assign(paste0(sample, "_tl"), tax_vec)
}
# Print the names of the created objects
ls(pattern = "_tl")

# Obtain the intersection of all
i_h2h3h4 <- Intersect(treatments_genera_list)
i_h2h3h4_d <- as.data.frame(i_h2h3h4, row.names = NULL, col.names = NULL) 
nrow(i_h2h3h4_d)
write.table(i_h2h3h4_d, "core_hydro_genera_otus_num.tsv", sep = "\t", row.names = FALSE)

```

#  Relaxed hydroponics core calculation


```R
# Convert lists of characters to numeric
hydro2_1_e_tl <- as.numeric(hydro2_1_e_tl)
hydro2_2_e_tl <- as.numeric(hydro2_2_e_tl)
hydro3_1_e_tl <- as.numeric(hydro3_1_e_tl)
hydro3_2_e_tl <- as.numeric(hydro3_2_e_tl)
hydro4_1_e_tl <- as.numeric(hydro4_1_e_tl)
hydro4_2_e_tl <- as.numeric(hydro4_2_e_tl)
hydro2_1_r_tl <- as.numeric(hydro2_1_r_tl)
hydro2_2_r_tl <- as.numeric(hydro2_2_r_tl)
hydro3_1_r_tl <- as.numeric(hydro3_1_r_tl)
hydro3_2_r_tl <- as.numeric(hydro3_2_r_tl)
hydro4_1_r_tl <- as.numeric(hydro4_1_r_tl)
hydro4_2_r_tl <- as.numeric(hydro4_2_r_tl)

# Create a list with the names of the lists
list_names <- c("hydro2_1_e_tl", "hydro2_2_e_tl", 
                "hydro3_1_e_tl", "hydro3_2_e_tl",
                "hydro4_1_e_tl", "hydro4_2_e_tl",
                "hydro2_1_r_tl", "hydro2_2_r_tl", 
                "hydro3_1_r_tl", "hydro3_2_r_tl",
                "hydro4_1_r_tl", "hydro4_2_r_tl")

# Convert list names into a list of lists
list_of_lists <- list(hydro2_1_e_tl, hydro2_2_e_tl, 
                hydro3_1_e_tl, hydro3_2_e_tl,
                hydro4_1_e_tl, hydro4_2_e_tl,
                hydro2_1_r_tl, hydro2_2_r_tl, 
                hydro3_1_r_tl, hydro3_2_r_tl,
                hydro4_1_r_tl, hydro4_2_r_tl)

# Function to get intersection of a list of lists
Intersect <- function(lists) {
  Reduce(intersect, lists)
}

# Loop to exclude a different list in each iteration
for (i in 1:length(list_of_lists)) {
  # Create a list excluding the list at position i
  temp_list <- list_of_lists[-i]
  # Calculate the intersection of the remaining lists
  intersection_result <- Intersect(temp_list)
  # Create the name for the output object
  output_name <- paste0("less_", list_names[i])
  # Save the result in the global environment with the created name
  assign(output_name, intersection_result)
}

# Print the names of the created objects
print(ls(pattern = "less_"))

list_core_relax <- c("less_hydro2_1_e_tl", "less_hydro2_1_r_tl", "less_hydro2_2_e_tl",
                    "less_hydro2_2_r_tl", "less_hydro3_1_e_tl", "less_hydro3_1_r_tl", 
                    "less_hydro3_2_e_tl", "less_hydro3_2_r_tl", "less_hydro4_1_e_tl", 
                    "less_hydro4_1_r_tl", "less_hydro4_2_e_tl", "less_hydro4_2_r_tl")

# Save relaxed core lists to files based on missing samples
for (i in list_core_relax) {
    file_name <- paste0(i, ".csv")
    lista <- get(i)
    dfl <- as.data.frame(lista)
    write.table(dfl, file_name, sep="\t", row.names = FALSE, col.names = FALSE)
}
```

# Plot core and relaxed core abundaces


```R
# Load core data
hydro_core <- as.vector(i_h2h3h4_d$i_h2h3h4)
hydro_core

# Load relaxed core data
hydro_core_relax <- read.table("hydro_core_relax.txt", header = FALSE)
hydro_core_relax <- as.vector(hydro_core_relax$V1)
hydro_core_relax

# hydro_core_relax.txt contains the genera of the relaxed core

# Combine core and relaxed core data
cores <- c(hydro_core, hydro_core_relax)
cores
length(cores)

# Filter samples to include only those with treatments
# Filtering to include only treatment samples
hydro_treat <- filter_taxa(
  subset_samples(tom_hydro_genus_ms, treatment_c == "treatment"),
  function(x) {sum(x > 0) > 0},
  prune = TRUE
)

# Subset taxa to include only those in the core and relaxed core lists
core_hydro_re <- subset_taxa(
  hydro_treat,
  rownames(tax_table(hydro_treat)) %in% cores
)
core_hydro_re

# Estimate relative abundance
# Transform counts to relative abundance
core_hydro_re_rel <- transform_sample_counts(
  core_hydro_re,
  function(x) x / sum(x)
)
core_hydro_re_rel
```


```R
#Obtain shared genera to determine the most abundant and then sort them in the heatmap
#extract data from the plot
p <- plot_bar(core_hydro_re_rel, "Genus")
yx <- p$data
yx <- as.data.frame(yx)

#sort data by abundance in descending order
yx <- yx[order(-yx$Abundance),]

#select desired columns
columnas_deseadas <- yx[, c("Abundance", "Genus")]

#list of target genera
core_treat <- c("_g__Gemmobacter", "_g__Kaistia", "_g__Caulobacter", "_g__Rhodobacter", 
                "_g__Luteolibacter", "_g__Pseudomonas", "_g__Dongia", "_g__Runella", 
                "_g__Leptothrix", "_g__Hyphomicrobium", "_g__Shinella", "_g__Devosia", 
                "_g__Aminobacter", "_g__Pseudoxanthomonas", "_g__Brevundimonas", "_g__", 
                "_g__Acinetobacter", "_g__Blastomonas", "_g__Dyadobacter", "_g__Sphingopyxis", 
                "_g__Novosphingobium", "_g__Flavobacterium", "_g__Caedibacter", "_g__Prosthecobacter", 
                "_g__Bosea", "_g__Hydrogenophaga", "_g__Ancylobacter", "_g__Neorhizobium", 
                "_g__Aurantimonas", "_g__Sphingomonas", "_g__Sphingobium",
                "_g__Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", 
                "_g__Defluviimonas", "_g__Bacillus")

#filter rows based on the genera
compartidos_table <- subset(columnas_deseadas, Genus %in% core_treat)

#calculate the sum of abundances by genus
suma_por_genero <- aggregate(Abundance ~ Genus, data = compartidos_table, sum)

#sort genera from lowest to highest sum
suma_por_genero_ordenada <- suma_por_genero[order(suma_por_genero$Abundance), ]
suma_por_genero_ordenada

#save the names of the genera in a vector
core_hydro_names <- as.vector(suma_por_genero_ordenada$Genus)
core_hydro_names
```


```R
#Obtain genera from the relaxed core to determine their abundance and then sort them in the heatmap
core_relaxed <- c("_g__Candidatus_Protochlamydia", "_g__Xanthobacter", "_g__Bradyrhizobium", "_g__Haematospirillum", 
                 "_g__Acidovorax", "_g__Defluviicoccus", "_g__Edaphobaculum", "_g__Nordella", "_g__Flavihumibacter",
                 "_g__Gemmatimonas", "_g__Hyphomonas", "_g__Porphyrobacter", "_g__Erythrobacter", 
                 "_g__Pleomorphomonas")

#filter rows based on the genera
table <- subset(columnas_deseadas, Genus %in% core_relaxed)

#calculate the sum of abundances by genus
suma_por_genero <- aggregate(Abundance ~ Genus, data = table, sum)

#sort genera from lowest to highest sum
suma_por_genero_ordenada <- suma_por_genero[order(suma_por_genero$Abundance), ]
suma_por_genero_ordenada

#save the names of the genera in a vector
core_relaxed_names <- as.vector(suma_por_genero_ordenada$Genus)
core_relaxed_names
```


```R
#heatmap core
#extract data from the plot
p <- plot_bar(core_hydro_re_rel, "Genus")
yx <- p$data
yx <- as.data.frame(yx)

#sort data by abundance in descending order
yx <- yx[order(-yx$Abundance),]

#order of the samples
order <- factor(yx$Sample, levels = c("hydro10", "hydro12", "hydro18", "hydro20", "hydro26", 
                                      "hydro28", "hydro9", "hydro11", "hydro17", "hydro19", 
                                      "hydro25", "hydro27"))

#order of the genera based on relaxed and core genera
order_genera_f <- factor(yx$Genus, levels = c(core_relaxed_names, core_hydro_names))

#plot relative abundance of phyla
ggplot(yx, aes(x = order, y = order_genera_f, fill = Abundance)) + 
  geom_raster() + 
  theme_minimal() + 
  theme(text = element_text(size = 20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient(low = "#bbd9fc", high = "#133b6b", na.value = "white", trans = "log2")
```

# Comparison of hydroponic and soil tomato genera


```R
#filtering rhizosphere 
rhizophere_hydro <- filter_taxa(subset_samples(treatment_genera, Rh_En == "rhizosphere"), function(x) {sum(x > 0) > 0}, prune = TRUE)
rhizophere_hydro

#filtering endosphere 
endosphere_hydro <- filter_taxa(subset_samples(treatment_genera, Rh_En == "endosphere"), function(x) {sum(x > 0) > 0}, prune = TRUE)
endosphere_hydro

#save rhizosphere hydroponics list
rhiz_hydro <- as.data.frame(tax_table(rhizophere_hydro))
rhiz_hydro_t <- apply(rhiz_hydro, 1, function(x) paste(x[-1], collapse = " | "))                                   
write.table(rhiz_hydro_t, "generos_rhizo_hydro.tsv", sep = "\t", row.names = FALSE)

#save endosphere hydroponics list
endo_hydro <- as.data.frame(tax_table(endosphere_hydro))
endo_hydro_t <- apply(endo_hydro, 1, function(x) paste(x[-1], collapse = " | "))                                   
write.table(endo_hydro_t, "generos_endo_hydro.tsv", sep = "\t", row.names = FALSE)

```


```R
#load genera data from hydroponic tomato
rhiz_hydro <- read.table("generos_rhizo_hydro.tsv", header = TRUE, sep = "\t")
endo_hydro <- read.table("generos_endo_hydro.tsv", header = TRUE, sep = "\t")

#load genera data from soil-grown tomato
rhiz_soil <- read.table("generos_rhizo_soil.tsv", header = TRUE, sep = "\t")
endo_soil <- read.table("generos_endo_soil.tsv", header = TRUE, sep = "\t")

#convert data to vectors
rhiz_hydro_v <- as.vector(rhiz_hydro$x)
endo_hydro_v <- as.vector(endo_hydro$x)
rhiz_soil_v <- as.vector(rhiz_soil$x)
endo_soil_v <- as.vector(endo_soil$x)

#create Venn diagram comparing rhizosphere and endosphere in hydroponics and soil
x2 <- list(rhiz_soil_v, rhiz_hydro_v, endo_soil_v, endo_hydro_v)

venn_hydro <- display_venn(x2, 
                           category.names = c("rhiz_soil_v", "rhiz_hydro_v", "endo_soil_v", "endo_hydro_v"),
                           lwd = 3, cex = 3,
                           lty = 'blank',
                           fill = c("#6f4040ff", "#3f8bc0ff", "#c8b6b6af", "#c0ddeba6"))
```

# Tomato core microbiome calculations


```R
#filtering treatment 
treatment_genera<-filter_taxa(subset_samples(tom_hydro_genus_ms, treatment_c=="treatment"), function (x) {sum(x > 0) > 0}, prune=TRUE)
treatment_genera

```


```R
# Construct an OTU table with genera to create an UpSet plot of tomatoes in hydroponics
tax_table_genus <- as.data.frame(tax_table(hydro_treat))

# Create a single line for the entire taxonomy by concatenating the taxonomy levels, excluding the first column
tt_genus_cat <- apply(tax_table_genus, 1, function(x) paste(x[-1], collapse = " | "))

# Convert the concatenated taxonomy into a data frame
tt_genus_cat_v <- as.data.frame(tt_genus_cat)

# Extract the OTU table from the hydro_treat object and convert it to a data frame
ott <- as.data.frame(otu_table(hydro_treat))

# Generate an OTU table with genera, but include the complete taxonomy line
# Merge the OTUs with their corresponding genera (full taxonomy line)
ot_genus <- cbind(ott, tt_genus_cat_v)

# Remove existing row names and set new row names based on the full taxonomy line
rownames(ot_genus) <- NULL
rownames(ot_genus) <- ot_genus$tt_genus_cat

# Remove the temporary full taxonomy column as it is now set as row names
ot_genus$tt_genus_cat <- NULL
```


```R
# Prepare the presence-absence table
entradaUpset <- ot_genus
# Convert all values greater than 0 to 1, indicating presence
entradaUpset[entradaUpset > 0] <- 1

# Save the presence-absence table to a file
write.table(entradaUpset, "tomatohydro_genera.tsv", sep = "\t", row.names = TRUE)

# Load the OTU table at the genus level for hydroponic tomato
genustable_hydrotom <- as.data.frame(read.table("tomatohydro_genera.tsv", header = TRUE, sep = "\t"))

# Load the OTU table at the genus level for soil-grown tomato
genustable_soiltom <- as.data.frame(read.table("tomatosoil_genera.tsv", header = TRUE, sep = "\t"))

# Add a new column 'RowNames' to store the row names
genustable_hydrotom$RowNames <- rownames(genustable_hydrotom)
genustable_soiltom$RowNames <- rownames(genustable_soiltom)

# Merge the tables based on the 'RowNames' column
entradaup_hydro_soil <- merge(genustable_hydrotom, genustable_soiltom, by = "RowNames", all = TRUE)

# Set the row names of the merged table based on the 'RowNames' column
rownames(entradaup_hydro_soil) <- entradaup_hydro_soil$RowNames

# Remove the 'RowNames' column as it's now set as the row names
entradaup_hydro_soil$RowNames <- NULL

# Replace NA values with 0 to indicate absence
entradaup_hydro_soil[is.na(entradaup_hydro_soil)] <- 0

# Define the sets for UpSet plot
all_hydro_soil <- c("JAL2RFT", "GTO2RFT", "JAL1RFT", "JAL3RFT", "JAL4RFT", "JAL5RFT",
                    "GTO1RFT", "SLP1RFT", "DGO1RFT", "ZAC1RFT", "NAY2RFT", "SIN2RFT",
                    "SIN1RFT", "AGS1RFT", "NAY3RFT", "GTO3RFT", "JAL2ECT", "SIN2ECT", 
                    "AGS1ECT", "GTO2ECT", "SLP1ECT", "ZAC1ECT", "DGO1ECT", "JAL3ECT", 
                    "NAY3ECT", "JAL1ECT", "JAL4ECT", "JAL5ECT", "GTO3ECT", "NAY2ECT", 
                    "SIN1ECT", "GTO1ECT", "hydro10", "hydro11", "hydro12", "hydro17",
                    "hydro18", "hydro19", "hydro20", "hydro25", 
                    "hydro26", "hydro27", "hydro28", "hydro9")

all_hydro <- c("hydro10", "hydro11", "hydro12", "hydro17",
               "hydro18", "hydro19", "hydro20", "hydro25", 
               "hydro26", "hydro27", "hydro28", "hydro9")

# Generate the UpSet plot to visualize intersections
tomateUpset <- upset(as.data.frame(entradaup_hydro_soil),
                     sets = all_hydro_soil, nintersects = 15, 
                     order.by = "freq", # Order by frequency (set size)
                     number.angles = 0, 
                     point.size = 2.5, 
                     line.size = 0.5, 
                     mainbar.y.label = "Genus intersections", 
                     sets.x.label = "Genus per sample",  
                     text.scale = c(4, 4, 4, 4, 3, 4)
)
tomateUpset
```

# Tomato core microbiome relative abundances


```R
# Generate the core tomato microbiome for both hydroponic and soil environments
core_hydro_soil_list <- c("_g__Prosthecobacter", "_g__Bacillus", "_g__Devosia", 
                          "_g__Sphingomonas", "_g__Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
                          "_g__Caulobacter", "_g__Dongia", "_g__Novosphingobium", "_g__Sphingobium")

# Subset the taxa to obtain only the core genera in hydroponics
core_hydro_soil <- subset_taxa(hydro_treat, Genus %in% core_hydro_soil_list)
core_hydro_soil

# Merge samples by treatment category
core_hydro_soil_m <- merge_samples(core_hydro_soil, "treatment_c")

# Calculate relative abundance for the merged samples
core_hydro_soil_m_r <- transform_sample_counts(core_hydro_soil_m, function(x) x / sum(x))
core_hydro_soil_m_r

# Plot the CORE heatmap in hydroponics 
# Extract data for the plot
p <- plot_bar(core_hydro_soil_m_r, "Genus")
yx <- p$data
yx <- as.data.frame(yx)

# Order the data by abundance in descending order
yx <- yx[order(-yx$Abundance),]

# Plot the relative abundance of core genera in hydroponics as a heatmap
ggplot(yx, aes(x = Sample, y = reorder(Genus, Abundance), fill = Abundance)) + 
  geom_raster() + 
  theme_minimal() + 
  theme(text = element_text(size = 18)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient(low = "#bbd9fc", high = "#133b6b", na.value = "white", trans = "log2")

```


```R
# Load the phyloseq object for tomato samples grown in soil
otu <- as.matrix(read.table("allsoil.otu", header = TRUE, row.names = 1))
# Create the OTU table object with taxa as rows
OTU <- otu_table(otu, taxa_are_rows = TRUE)
# Load the taxonomy table from the file and convert it to a matrix
taximat <- as.matrix(read.table("tax_soil.clean", header = TRUE, row.names = 1))
# Create the taxonomy table object
taxi <- tax_table(taximat)
# Combine the OTU and taxonomy tables into a phyloseq object for soil-grown tomatoes
tomate_soil <- phyloseq(OTU, taxi)
# Load the sample metadata from the file
sample_data <- read.table("metadata2_o.tsv", header = TRUE, row.names = 1, sep = "\t")
# Create a sample data object with the necessary columns from the metadata
sampledata <- sample_data(data.frame(kind = sample_data$kind, 
                                     tag = sample_data$tag, 
                                     type = sample_data$type, 
                                     row.names = sample_names(tomate_soil)))
barajas2020 <- phyloseq(OTU, sampledata, taxi)
# Display the final phyloseq object
barajas2020
```


```R
# Perform taxonomic agglomeration at the genus level
tomate_soil_glom <- tax_glom(barajas2020, taxrank = "genus")
tomate_soil_glom
```


```R
# Subset the taxa to obtain only the core genera in soil
core_soil <- subset_taxa(tom_soil, genus %in% core_hydro_soil_list)
core_soil

# Merge the samples by type
core_soil_m <- merge_samples(core_soil, "type")

# Estimate relative abundance for the merged samples
core_soil_r <- transform_sample_counts(core_soil_m, function(x) x / sum(x))
core_soil_r

# Extract data for plotting the CORE heatmap
p <- plot_bar(core_soil_r, "genus")
yx <- p$data
yx <- as.data.frame(yx)

# Order the data by abundance in descending order
yx <- yx[order(-yx$Abundance),]

# Plot the relative abundance of core genera in soil as a heatmap
ggplot(yx, aes(x = Sample, y = reorder(genus, Abundance), fill = Abundance)) + 
  geom_raster() + 
  theme_minimal() + 
  theme(text = element_text(size = 18)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_gradient(low = "#ccc3b8", high = "#6f4040bf", na.value = "white", trans = "log2")
```
