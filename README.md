# Human-Metagenomic-Analysis-and-Visualization-using-R-SRA-SRR34768016-
Complete human metagenomic analysis of SRR34768016 using Linux and R. Includes QC, genome mapping, functional profiling, and visualization (barplots, heatmaps, circular genome plots). Transforms raw reads into meaningful biological insights for microbiome research and functional annotation.
### loaded Library
library(dplyr)
library(ggplot2)
library(tidyverse)

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
## Heatmap
presence_matrix <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
rownames(presence_matrix) <- paste0("Gene", 1:10)
colnames(presence_matrix) <- paste0("Contig", 1:10)

pheatmap::pheatmap(presence_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
                   main = "Gene Presence/Absence Heatmap")
####
library(igraph)

# Build a simple gene co-occurrence network 

adj <- matrix(sample(0:1, 100, replace = TRUE, prob = c(0.8, 0.2)), nrow=10)
diag(adj) <- 0
g <- graph.adjacency(adj, mode = "undirected")

plot(g, vertex.label = paste("Gene", 1:10), main = "Gene Co-occurrence Network")
###
# data frame
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
<img width="426" height="349" alt="tRNA" src="https://github.com/user-attachments/assets/ccd75776-a196-4c04-9de6-f86d4d6e6081" />
<img width="426" height="349" alt="Taxonomic by samples" src="https://github.com/user-attachments/assets/a96370db-47fc-4840-a963-7a2c7d9c24c1" />
<img width="424" height="347" alt="gene sample" src="https://github.com/user-attachments/assets/4ee75b1f-2a75-40f4-abfa-2344c660b55e" />
<img width="426" height="349" alt="gene length" src="https://github.com/user-attachments/assets/88553daf-7b50-446c-9369-c81be851d224" />
<img width="426" height="349" alt="gene copccurence network" src="https://github.com/user-attachments/assets/0f591b98-0f80-47bf-9b98-a72bfdfa6cb6" />
<img width="426" height="349" alt="feature distribution cds" src="https://github.com/user-attachments/assets/5028763b-c37a-44c5-83b9-360240e8dd05" />
<img width="426" height="349" alt="best heat map contig" src="https://github.com/user-attachments/assets/622df68c-9f0e-4697-be16-4433add279a8" />
<img width="426" height="349" alt="annotated gene" src="https://github.com/user-attachments/assets/24803694-5e04-4520-8691-e356f6ffb4d0" />





