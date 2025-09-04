# Load necessary libraries
library(dplyr)
library(ggplot2)

######last updated
library(tidyverse)

# Set your working directory where sample1.tsv is located
setwd("~/human meta")

# Load the Prokka annotation table
prokka_data <- read_tsv("sample1.tsv")

# Summarize counts of each gene function (product)
gene_counts <- prokka_data %>%
  count(product, sort = TRUE)

cat("Top gene functions (including hypothetical):\n")
print(gene_counts)

# Filter out hypothetical proteins to focus on known gene functions
known_genes <- gene_counts %>%
  filter(!str_detect(product, regex("hypothetical", ignore_case = TRUE)))

cat("\nKnown gene functions:\n")
print(known_genes)

# Plot the known gene functions if any exist
if (nrow(known_genes) > 0) {
  ggplot(known_genes, aes(x = reorder(product, n), y = n)) +
    geom_col(fill = "darkgreen") +
    coord_flip() +
    labs(title = "Known Gene Functions",
         x = "Gene Function",
         y = "Count") +
    theme_minimal()
} else {
  cat("\nNo known gene functions to plot (mostly hypothetical proteins).\n")
}

# Summarize feature types using the 'ftype' column from your data
feature_counts <- prokka_data %>%
  count(ftype, sort = TRUE)

cat("\nFeature type counts:\n")
print(feature_counts)

# Plot feature type distribution
ggplot(feature_counts, aes(x = reorder(ftype, n), y = n)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(title = "Feature Type Distribution",
       x = "Feature Type",
       y = "Count") +
  theme_minimal()

# Save the gene and feature counts as CSV for reference
write_csv(gene_counts, "gene_function_counts.csv")
write_csv(feature_counts, "feature_type_counts.csv")

cat("\nSummary tables saved as 'gene_function_counts.csv' and 'feature_type_counts.csv'\n")

##
# For org.Hs.eg.db (human annotation package)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c( "org.Hs.eg.db"))

# For pheatmap from CRAN
install.packages("pheatmap")

# Load necessary libraries
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)  # change or remove based on your organism

# Set working directory
setwd("~/human meta")

# Load Prokka annotation table
prokka_data <- read_tsv("sample1.tsv")

# 1) Summarize gene function counts (product)
gene_counts <- prokka_data %>%
  count(product, sort = TRUE)

# Filter out hypothetical proteins for meaningful analysis
known_genes <- prokka_data %>%
  filter(!str_detect(product, regex("hypothetical", ignore_case = TRUE)))

# 2) Functional category summary (assuming COG available)
if ("COG" %in% colnames(prokka_data)) {
  cog_counts <- prokka_data %>%
    filter(!is.na(COG)) %>%
    count(COG, sort = TRUE)
  
  print("Top COG categories:")
  print(head(cog_counts, 10))
} else {
  cat("No COG column found; skipping COG enrichment.\n")
}

# 3) Prepare for enrichment analysis (example with GO terms, adjust for your data)

# If you have gene IDs and GO annotations, you can do enrichment with clusterProfiler
# Here we fake a gene list for demo - replace this with real gene IDs and categories

# Example: Convert gene products to Entrez IDs (if possible)
# entrez_ids <- mapIds(org.Hs.eg.db, keys=known_genes$product, column="ENTREZID", keytype="SYMBOL", multiVals="first")
# entrez_ids <- na.omit(entrez_ids)

# Run GO enrichment - example for Biological Process
# ego <- enrichGO(gene = entrez_ids,
#                 OrgDb = org.Hs.eg.db,
#                 ont = "BP",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff = 0.05,
#                 readable = TRUE)
# print(head(ego))

# 4) Visualization: Heatmap of top gene functions

# Create matrix for heatmap: products vs counts
top_products <- gene_counts %>% head(20)
heatmap_data <- matrix(top_products$n, nrow = 1)
rownames(heatmap_data) <- "Sample1"
colnames(heatmap_data) <- top_products$product

pheatmap(heatmap_data,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         main = "Top 20 Gene Functions Heatmap",
         angle_col = 45,
         fontsize_col = 10)

# 5) Export gene list for pathway analysis

# Export known genes to CSV
write_csv(known_genes, "known_genes_for_pathway_analysis.csv")
cat("Exported known genes for pathway analysis.\n")

####----------
library(ggplot2)

ggplot(prokka_data, aes(x = length_bp)) +
  geom_histogram(binwidth = 50, fill = "purple", color = "black") +
  labs(title = "Distribution of Gene Lengths",
       x = "Gene Length (bp)",
       y = "Count") +
  theme_minimal()
##
if("ftype" %in% colnames(prokka_data)) {
  feature_counts <- prokka_data %>% count(ftype, sort = TRUE)
  ggplot(feature_counts, aes(x = reorder(ftype, n), y = n)) +
    geom_col(fill = "green") +
    coord_flip() +
    labs(title = "Feature Type Distribution",
         x = "Feature Type",
         y = "Count") +
    theme_minimal()
} else {
  message("No 'ftype' column found.")
}
##
library(Biostrings)
library(ggplot2)

# Load contigs fasta file
contigs <- readDNAStringSet("contigs.fasta")

# Get contig lengths
contig_lengths <- width(contigs)

# Summary stats
summary(contig_lengths)

# Plot length distribution
ggplot(data.frame(length = contig_lengths), aes(x = length)) +
  geom_histogram(binwidth = 500, fill = "blue", color = "black") +
  scale_x_log10() +
  labs(title = "Contig Length Distribution", x = "Contig Length (log scale)", y = "Count") +
  theme_minimal()
####
# Load required packages
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("Biostrings")
}
library(Biostrings)
library(ggplot2)

# Read contigs fasta file
contigs <- readDNAStringSet("contigs.fasta")

# Calculate contig lengths
contig_lengths <- width(contigs)

# Summary statistics of contig lengths
print(summary(contig_lengths))

# Plot contig length distribution
ggplot(data.frame(length = contig_lengths), aes(x = length)) +
  geom_histogram(binwidth = 500, fill = "steelblue", color = "black") +
  scale_x_log10() +  # Log scale for better visualization
  labs(title = "Contig Length Distribution",
       x = "Contig Length (log scale)",
       y = "Count") +
  theme_minimal()
###
library(Biostrings)

# Calculate GC content
gc_content <- letterFrequency(contigs, letters = c("G", "C"), as.prob = TRUE)
gc_percent <- rowSums(gc_content) * 100

# Plot GC content distribution
ggplot(data.frame(gc = gc_percent), aes(x = gc)) +
  geom_histogram(binwidth = 1, fill = "yellow", color = "black") +
  labs(title = "GC Content Distribution per Contig",
       x = "GC Content (%)",
       y = "Count") +
  theme_minimal()
###
df <- data.frame(length = contig_lengths, gc = gc_percent)

ggplot(df, aes(x = length, y = gc)) +
  geom_point(alpha = 0.5, color = "orange") +
  scale_x_log10() +
  labs(title = "Contig Length vs GC Content",
       x = "Contig Length (log scale)",
       y = "GC Content (%)") +
  theme_minimal()
###
library(dplyr)

# Assuming prokka_data is loaded with a column ftype representing feature types
feature_counts <- prokka_data %>%
  count(ftype)

ggplot(feature_counts, aes(x = "", y = n, fill = ftype)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "Feature Type Composition") +
  theme_void() +
  theme(legend.title = element_blank())
###
top_genes <- prokka_data %>%
  filter(!grepl("hypothetical", product, ignore.case = TRUE)) %>%
  count(product, sort = TRUE) %>%
  head(10)

ggplot(top_genes, aes(x = reorder(product, n), y = n)) +
  geom_col(fill = "coral") +
  coord_flip() +
  labs(title = "Top 10 Annotated Genes",
       x = "Gene Product",
       y = "Count") +
  theme_minimal()
##
# Example dummy matrix (replace with real presence/absence data)
presence_matrix <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
rownames(presence_matrix) <- paste0("Gene", 1:10)
colnames(presence_matrix) <- paste0("Contig", 1:10)

pheatmap::pheatmap(presence_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
                   main = "Gene Presence/Absence Heatmap")
####
library(igraph)

# Build a simple gene co-occurrence network example
# (You'd create an adjacency matrix from your real data)
adj <- matrix(sample(0:1, 100, replace = TRUE, prob = c(0.8, 0.2)), nrow=10)
diag(adj) <- 0
g <- graph.adjacency(adj, mode = "undirected")

plot(g, vertex.label = paste("Gene", 1:10), main = "Gene Co-occurrence Network")
###
# Example data frame
taxa <- data.frame(
  Sample = rep(c("Sample1", "Sample2"), each = 4),
  Taxon = rep(c("Firmicutes", "Proteobacteria", "Bacteroidetes", "Actinobacteria"), 2),
  Abundance = c(30, 50, 10, 10, 40, 30, 20, 10)
)

ggplot(taxa, aes(x = Sample, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity") +
  labs(title = "Taxonomic Composition by Sample", y = "Relative Abundance (%)") +
  theme_minimal()
###
install.packages("sunburstR")
install.packages("d3r")

library(d3r)
library(sunburstR)

# Example hierarchical data (replace with your real functional categories)
df <- data.frame(
  path = c("Metabolism|Carbohydrate metabolism|Glycolysis",
           "Metabolism|Amino acid metabolism|Glutamate",
           "Information processing|Replication|DNA polymerase"),
  value = c(50, 30, 20)
)

sb_data <- d3_nest(df, value_cols = "value", root = "Root")
sunburst(sb_data)
###
install.packages("plotly")  # if not installed
library(plotly)

mat <- matrix(rnorm(100), nrow = 10)
rownames(mat) <- paste("Gene", 1:10)
colnames(mat) <- paste("Sample", 1:10)

plot_ly(
  z = mat,
  x = colnames(mat),
  y = rownames(mat),
  type = "heatmap",
  colorscale = "Viridis"
) %>%
  layout(title = "Interactive Heatmap")


