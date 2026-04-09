
################################################################################
# Microbiome Analysis Demo: Diversity, Networks, and Differential Abundance
################################################################################

# --- Step 0: Install and Load Libraries ---
# Installs and loads all required packages for microbiome analysis
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("phyloseq")) install.packages("phyloseq")
if (!require("vegan")) install.packages("vegan")
if (!require("ggpubr")) install.packages("ggpubr")
if (!require("picante")) install.packages("picante")
if (!require("svglite")) install.packages("svglite")
if (!require("ALDEx2")) BiocManager::install("ALDEx2")

library(tidyverse)
library(phyloseq)
library(vegan)
library(ggpubr)
library(picante)
library(svglite)
library(ALDEx2)

# --- Step 1: Data Import & Cleaning ---
# Reads OTU, taxonomy, metadata, and phylogenetic tree, then merges into a phyloseq object
otu_table_raw <- read.csv("otu_table.csv", row.names = 1) %>% as.matrix()
taxa_table_raw <- read.csv("taxa_table.csv", row.names = 1) %>% as.matrix()
metadata_raw <- read.csv("metadata.csv", sep = '\t', row.names = 1)
agp_dc_tree <- read_tree("tree.nwk")

physeq <- phyloseq(
  otu_table(otu_table_raw, taxa_are_rows = TRUE),
  tax_table(taxa_table_raw),
  sample_data(metadata_raw),
  phy_tree(agp_dc_tree)
)

# Filters samples with >1000 reads
physeq_1k <- prune_samples(sample_sums(physeq) > 1000, physeq)
print(physeq_1k)

# --- Step 2: Alpha Diversity ---
# Calculates richness and diversity indices, then plots boxplots comparing groups
alpha_div <- estimate_richness(physeq_1k, measures = c("Observed", "Shannon", "Simpson"))
plot_richness(physeq_1k, x="Subset", color = "Subset", 
              measures = c("Shannon", "Observed", "Simpson"), 
              title = "Combined Alpha Diversity measures") + 
  geom_boxplot() + theme_bw() + theme(axis.text = element_text(size = 6))

# Performs Wilcoxon test between Diabetic vs Healthy for observed richness
comparisons <- list(c("Diabetic", "Healthy"))
plot_richness(physeq_1k, x="Subset", color = "Subset", measures = "Observed", title = "Observed Richness") +
  geom_boxplot() + stat_compare_means(comparisons = comparisons, method = "wilcox.test") + 
  theme_bw() + theme(axis.text = element_text(size = 6))

# --- Step 3: Beta Diversity ---
# Computes UniFrac distances, performs PCoA ordination, plots clustering with ellipses
ord_unifrac <- ordinate(physeq_1k, method = "PCoA", distance = "unifrac")
plot_ordination(physeq_1k, ord_unifrac, color = "Subset", title = "Unifrac PCoA") +
  stat_ellipse() + theme_bw() +
  theme(legend.key.size = unit(0.1,'cm'), legend.text = element_text(size = 7))

# Tests variance homogeneity (beta dispersion) and runs PERMANOVA
beta_dist <- phyloseq::distance(physeq_1k, method = "unifrac")
beta_disper_res <- betadisper(d = beta_dist, group = sample_data(physeq_1k)$Subset)
anova(beta_disper_res)

adonis2(beta_dist ~ Subset + age_cat, data = data.frame(sample_data(physeq_1k)), 
        permutations = 999, by = "terms")

# --- Step 4: Canonical Correspondence Analysis (CCA) ---
# Constrained ordination to see how metadata explains OTU distribution
cca_res <- ordinate(physeq_1k, method = "CCA", formula = ~ Subset + age_cat)
cca_signif <- anova.cca(cca_res, permutations = 999)
print(cca_signif)

# --- Step 5: Differential Abundance (ALDEx2) ---
# Runs ALDEx2 to test differential abundance between groups
counts_matrix <- as.matrix(otu_table(physeq_1k))
group_vector <- as.character(sample_data(physeq_1k)$Subset)
if (ncol(counts_matrix) != length(group_vector)) {
  counts_matrix <- t(counts_matrix)
}
x <- aldex(counts_matrix, group_vector, mc.samples = 128, test = "t", effect = TRUE, denom = "all")

# Volcano plot of effect size vs p-value, highlights significant taxa
plot(x$effect, -log10(x$we.ep), pch = 19, cex = 0.5, col = "darkgrey",
     xlab = "Effect Size", ylab = "-log10(p-value)",
     main = "ALDEx2: Diabetic vs Healthy")
sig_taxa <- x[x$we.eBH < 0.05, ]
points(sig_taxa$effect, -log10(sig_taxa$we.ep), col = "red", pch = 19, cex = 0.7)
abline(h = -log10(0.05), lty = 2, col = "blue")
abline(v = 0, lty = 1, col = "black")

# --- Step 6: Network Analysis ---
# Builds a co-occurrence network of taxa for Diabetic samples
physeq_species <- tax_glom(physeq_1k, taxrank = "Species")
physeq_sp_flt <- filter_taxa(physeq_species, function(x) sum(x) > 10, TRUE)
tax <- data.frame(tax_table(physeq_sp_flt))
taxa_names(physeq_sp_flt) <- tax$Species

diabetes_ps <- phyloseq::subset_samples(physeq_sp_flt, Subset %in% c('Diabetic'))
