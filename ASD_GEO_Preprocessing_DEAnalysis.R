# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GEOquery", "dplyr", "tibble"))
library(GEOquery)
library(dplyr)
library(tibble)

# Function to download and preprocess a GEO dataset by its GSE accession
preprocess_GEO_dataset <- function(gse_id) {
  message("Processing dataset: ", gse_id)
  # Download the GEO dataset (ExpressionSet)
  gse <- getGEO(gse_id, GSEMatrix = TRUE)
  # If multiple ExpressionSets, select the first
  gse <- gse[[1]]
  
  # Get platform ID and download platform annotation (GPL)
  platform_id <- annotation(gse)
  gpl <- getGEO(platform_id, AnnotGPL = TRUE)
  
  # Extract annotation table from GPL object
  gpl_annot <- Table(gpl)
  
  # Extract expression matrix and convert to data frame
  exprs_df <- exprs(gse) %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  
  # Select relevant columns from annotation (probe ID and gene symbol)
  # Make column names lowercase for easier matching
  colnames(gpl_annot) <- tolower(colnames(gpl_annot))
  
  # Try to find the column that contains gene symbols
  possible_gene_cols <- c("gene symbol", "gene_symbol", "symbol", "gene_symbols")
  gene_symbol_col <- possible_gene_cols[possible_gene_cols %in% colnames(gpl_annot)][1]
  
  if (is.na(gene_symbol_col)) {
    stop("Gene symbol column not found in platform annotation.")
  }
  
  gpl_annot_small <- gpl_annot %>%
    select(id, !!sym(gene_symbol_col)) %>%
    rename(ID = id, Gene.symbol = !!sym(gene_symbol_col))
  
  # Join expression data with gene symbols by probe ID
  exprs_annotated <- left_join(exprs_df, gpl_annot_small, by = "ID")
  
  # Filter out probes without gene symbols
  exprs_filtered <- exprs_annotated %>%
    filter(!is.na(Gene.symbol), Gene.symbol != "")
  
  message("Finished processing ", gse_id, 
          ": ", nrow(exprs_filtered), " probes with gene symbols retained.")
  
  return(exprs_filtered)
}

# Process each dataset and save results
datasets <- c("GSE26415", "GSE42133", "GSE111175")

for (ds in datasets) {
  processed_data <- preprocess_GEO_dataset(ds)
  
  # Save to CSV for later merging and batch correction
  write.csv(processed_data, paste0("exprs_", ds, "_annotated.csv"), row.names = FALSE)
}

exprs_26415 <- read.csv("exprs_GSE26415_annotated.csv")
exprs_42133 <- read.csv("exprs_GSE42133_annotated.csv")
exprs_111175 <- read.csv("exprs_GSE111175_annotated.csv")
str(exprs_26415)
str(exprs_42133)
str(exprs_111175)
# Check number of unique genes
length(unique(exprs_26415$Gene.symbol))
length(unique(exprs_42133$Gene.symbol))
length(unique(exprs_111175$Gene.symbol))

# Collapse multiple probes per gene by averaging expression values to get one row per gene per sample
collapse_genes <- function(exprs_data) {
  exprs_data %>%
    select(-ID) %>%
    group_by(Gene.symbol) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%
    ungroup()
}

exprs_26415_collapsed <- collapse_genes(exprs_26415)
exprs_42133_collapsed <- collapse_genes(exprs_42133)
exprs_111175_collapsed <- collapse_genes(exprs_111175)

# Keep only genes common to all datasets for consistent comparison and integration
common_genes <- Reduce(intersect, list(
  exprs_26415_collapsed$Gene.symbol,
  exprs_42133_collapsed$Gene.symbol,
  exprs_111175_collapsed$Gene.symbol
))

# Subset to common genes
exprs_26415_final <- exprs_26415_collapsed %>% filter(Gene.symbol %in% common_genes)
exprs_42133_final <- exprs_42133_collapsed %>% filter(Gene.symbol %in% common_genes)
exprs_111175_final <- exprs_111175_collapsed %>% filter(Gene.symbol %in% common_genes)

# Check dimensions
dim(exprs_26415_final)
dim(exprs_42133_final)
dim(exprs_111175_final)


# Remove Gene.symbol temporarily and transpose datasets
exprs_26415_mat <- as.data.frame(t(exprs_26415_final[,-1]))
colnames(exprs_26415_mat) <- exprs_26415_final$Gene.symbol

exprs_42133_mat <- as.data.frame(t(exprs_42133_final[,-1]))
colnames(exprs_42133_mat) <- exprs_42133_final$Gene.symbol

exprs_111175_mat <- as.data.frame(t(exprs_111175_final[,-1]))
colnames(exprs_111175_mat) <- exprs_111175_final$Gene.symbol

# Add a batch column
exprs_26415_mat$batch <- "GSE26415"
exprs_42133_mat$batch <- "GSE42133"
exprs_111175_mat$batch <- "GSE111175"

# Combine them
combined_data <- bind_rows(exprs_26415_mat, exprs_42133_mat, exprs_111175_mat)
write.csv(combined_data, "combined_expression_data.csv", row.names = TRUE)

#Load the GSE objects:
gse26415 <- getGEO("GSE26415", GSEMatrix = TRUE)[[1]]
options(timeout = 600)  # Increase timeout to 10 minutes
gse42133 <- getGEO("GSE42133", GSEMatrix = TRUE)[[1]]
gse111175 <- getGEO("GSE111175", GSEMatrix = TRUE)[[1]]

# Extract the sample metadata from each dataset
pheno_26415 <- pData(gse26415)
pheno_42133 <- pData(gse42133)
pheno_111175 <- pData(gse111175)
# Look at the column names
colnames(pheno_26415)
colnames(pheno_42133)
colnames(pheno_111175)

# Look at first few rows to guess where diagnosis is stored
head(pheno_26415)
head(pheno_42133)
head(pheno_111175)

# Get phenotype (sample metadata) gse42133
pheno <- pData(gse42133)

# Extract the diagnosis column
diagnosis <- pheno$`dx (diagnosis):ch1`

# Check values
table(diagnosis)

# Load the phenotype data gse26415
pheno_26415 <- pData(gse26415)

# Check unique labels
unique(pheno_26415$`sample type:ch1`)

# Load phenotype data
pheno_26415 <- pData(gse26415)

# Extract diagnosis column
diagnosis <- pheno_26415$`sample type:ch1`

# Keep only "ASD" and "control" samples
keep_idx <- diagnosis %in% c("ASD", "control")

# Filter phenotype data
pheno_26415_clean <- pheno_26415[keep_idx, ]

# Create clean labels: capitalize 'control'
labels_26415 <- diagnosis[keep_idx]
labels_26415[labels_26415 == "control"] <- "Control"

# add labels as a new column
pheno_26415_clean$CleanLabel <- labels_26415
# Making sure expression matrix columns are in same order
exprs_26415_filtered <- exprs_26415_final[, rownames(pheno_26415_clean)]
colnames(pheno_111175)
# Load phenotype data gse111175
pheno_111175 <- pData(gse111175)

# Identify the correct diagnosis column
diagnosis_111175 <- pheno_111175$`diagnosis:ch1`

# Filter for ASD and TD (consider TD as Control)
keep_111175 <- diagnosis_111175 %in% c("ASD", "TD")
pheno_111175_clean <- pheno_111175[keep_111175, ]

# Clean labels
labels_111175 <- diagnosis_111175[keep_111175]
labels_111175[labels_111175 == "TD"] <- "Control"
pheno_111175_clean$CleanLabel <- labels_111175

# Filter expression matrix
exprs_111175_filtered <- exprs_111175_final[, rownames(pheno_111175_clean)]

# Use the already clean phenotype data
pheno_42133_clean <- pData(gse42133)

# Extract labels
labels_42133 <- pheno_42133_clean$`dx (diagnosis):ch1`

# Ensure all "control" is written as "Control" 
labels_42133[labels_42133 == "control"] <- "Control"

# store in a column for consistency
pheno_42133_clean$CleanLabel <- labels_42133

# Subset expression matrix to matching samples
exprs_42133_filtered <- exprs_42133_final[, rownames(pheno_42133_clean)]


dim(exprs_26415_filtered)
dim(exprs_42133_filtered)
dim(exprs_111175_filtered)


combined_filtered_exprs <- bind_rows(
  as.data.frame(t(exprs_26415_filtered)),
  as.data.frame(t(exprs_42133_filtered)),
  as.data.frame(t(exprs_111175_filtered))
)
combined_labels <- c(labels_26415, labels_42133, labels_111175)


combined_metadata <- data.frame(
  SampleID = rownames(combined_filtered_exprs),
  Label = combined_labels,
  Batch = rep(c("GSE26415", "GSE42133", "GSE111175"),
              times = c(length(labels_26415), length(labels_42133), length(labels_111175)))
)



# Reattach sample IDs as row names for each filtered expression matrix
exprs_26415_filtered <- exprs_26415_final[, rownames(pheno_26415_clean)]
exprs_42133_filtered <- exprs_42133_final[, rownames(pheno_42133)]  # already clean
exprs_111175_filtered <- exprs_111175_final[, rownames(pheno_111175_clean)]

# Transpose and bind rows
transpose_and_annotate <- function(expr_mat, pheno, batch_id) {
  df <- as.data.frame(t(expr_mat))
  df$SampleID <- rownames(df)
  df$Label <- pheno$CleanLabel
  df$Batch <- batch_id
  return(df)
}

df_26415 <- transpose_and_annotate(exprs_26415_filtered, pheno_26415_clean, "GSE26415")
df_42133 <- transpose_and_annotate(exprs_42133_filtered, pheno_42133, "GSE42133")
df_111175 <- transpose_and_annotate(exprs_111175_filtered, pheno_111175_clean, "GSE111175")

# Combine all
combined_filtered_data <- bind_rows(df_26415, df_42133, df_111175)

# Save to CSV
write.csv(combined_filtered_data, "combined_ASD_Control_expression_data.csv", row.names = FALSE)


# Extract just metadata from each dataset
meta_26415 <- data.frame(
  CleanLabel = pheno_26415_clean$CleanLabel,
  SampleID = rownames(pheno_26415_clean),
  Batch = "GSE26415"
)


meta_42133 <- data.frame(
  CleanLabel = pheno_42133$`dx (diagnosis):ch1`,
  SampleID = rownames(pheno_42133),
  Batch = "GSE42133"
)

meta_111175 <- data.frame(
  CleanLabel = pheno_111175_clean$CleanLabel,
  SampleID = rownames(pheno_111175_clean),
  Batch = "GSE111175"
)

# Combine metadata
metadata_combined <- bind_rows(meta_26415, meta_42133, meta_111175)
write.csv(metadata_combined, "combined_sample_metadata.csv", row.names = FALSE)

# Flip expression matrix so each row is a sample and attach sample names.
exprs_26415_mat <- t(exprs_26415_final[,-1])
colnames(exprs_26415_mat) <- exprs_26415_final$Gene.symbol
exprs_26415_df <- as.data.frame(exprs_26415_mat)
exprs_26415_df$SampleID <- rownames(exprs_26415_df)

library(dplyr)

# Build metadata for each dataset

# GSE26415
pheno_26415 <- pData(gse26415)
labels_26415 <- pheno_26415$`sample type:ch1`
keep_26415 <- labels_26415 %in% c("ASD", "control")
labels_26415[labels_26415 == "control"] <- "Control"

meta_26415 <- data.frame(
  SampleID = rownames(pheno_26415)[keep_26415],
  Label = labels_26415[keep_26415],
  Batch = "GSE26415",
  stringsAsFactors = FALSE
)

# GSE42133
pheno_42133 <- pData(gse42133)
labels_42133 <- pheno_42133$`dx (diagnosis):ch1`
keep_42133 <- labels_42133 %in% c("ASD", "Control")  # already proper labels

meta_42133 <- data.frame(
  SampleID = rownames(pheno_42133)[keep_42133],
  Label = labels_42133[keep_42133],
  Batch = "GSE42133",
  stringsAsFactors = FALSE
)

# GSE111175
pheno_111175 <- pData(gse111175)
labels_111175 <- pheno_111175$`subject status:ch1`
keep_111175 <- labels_111175 %in% c("ASD", "TD")
labels_111175[labels_111175 == "TD"] <- "Control"

meta_111175 <- data.frame(
  SampleID = rownames(pheno_111175)[keep_111175],
  Label = labels_111175[keep_111175],
  Batch = "GSE111175",
  stringsAsFactors = FALSE
)
colnames(pheno_111175_clean)
unique(pheno_111175_clean$"dianosis:ch1")

# Filter expression matrices and transpose

filter_and_transpose <- function(exprs_final, meta_df) {
  filtered <- exprs_final[, c("Gene.symbol", meta_df$SampleID)]
  filtered <- filtered[!duplicated(filtered$Gene.symbol), ]
  gene_symbols <- filtered$Gene.symbol
  exprs_mat <- t(filtered[, -1])
  colnames(exprs_mat) <- gene_symbols
  exprs_mat <- as.data.frame(exprs_mat)
  exprs_mat$SampleID <- rownames(exprs_mat)
  return(exprs_mat)
}

exprs_26415_clean <- filter_and_transpose(exprs_26415_final, meta_26415)
exprs_42133_clean <- filter_and_transpose(exprs_42133_final, meta_42133)
exprs_111175_clean <- filter_and_transpose(exprs_111175_final, meta_111175)

# Join expression data with metadata

merged_26415 <- inner_join(meta_26415, exprs_26415_clean, by = "SampleID")
merged_42133 <- inner_join(meta_42133, exprs_42133_clean, by = "SampleID")
merged_111175 <- inner_join(meta_111175, exprs_111175_clean, by = "SampleID")

# Combine all datasets

combined_exprs <- bind_rows(merged_26415, merged_42133, merged_111175)

# Save to CSV

write.csv(combined_exprs,"combined_ASD_Control_expression.csv", row.names = FALSE)

exprs_26415 <- read.csv("C:\\Users\\GARIMA\\OneDrive\\Desktop\\dissertation\\ASD_Biomarker_discovery\\exprs_GSE26415_annotated.csv")
# Load annotated expression data using double backslashes in the path

exprs_42133 <- read.csv("C:\\Users\\GARIMA\\OneDrive\\Desktop\\dissertation\\ASD_Biomarker_discovery\\exprs_GSE42133_annotated.csv")
exprs_111175 <- read.csv("C:\\Users\\GARIMA\\OneDrive\\Desktop\\dissertation\\ASD_Biomarker_discovery\\exprs_GSE111175_annotated.csv")
library(dplyr)

collapse_genes <- function(df) {
  df %>%
    group_by(Gene.symbol) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
    ungroup()
}

exprs_26415_collapsed <- collapse_genes(exprs_26415)
exprs_42133_collapsed <- collapse_genes(exprs_42133)
exprs_111175_collapsed <- collapse_genes(exprs_111175)

common_genes <- Reduce(intersect, list(
  exprs_26415_collapsed$Gene.symbol,
  exprs_42133_collapsed$Gene.symbol,
  exprs_111175_collapsed$Gene.symbol
))

exprs_26415_filtered <- exprs_26415_collapsed %>% filter(Gene.symbol %in% common_genes)
exprs_42133_filtered <- exprs_42133_collapsed %>% filter(Gene.symbol %in% common_genes)
exprs_111175_filtered <- exprs_111175_collapsed %>% filter(Gene.symbol %in% common_genes)

transpose_exprs <- function(df, batch_name) {
  mat <- as.data.frame(t(df[,-1]))
  colnames(mat) <- df$Gene.symbol
  mat$SampleID <- rownames(mat)
  mat$Batch <- batch_name
  return(mat)
}

exprs_26415_mat <- transpose_exprs(exprs_26415_filtered, "GSE26415")
exprs_42133_mat <- transpose_exprs(exprs_42133_filtered, "GSE42133")
exprs_111175_mat <- transpose_exprs(exprs_111175_filtered, "GSE111175")
combined_exprs <- bind_rows(exprs_26415_mat, exprs_42133_mat, exprs_111175_mat)
write.csv(combined_exprs, "combined_expression_matrix.csv", row.names = FALSE)

library(GEOquery)
gse26415 <- getGEO("GSE26415", GSEMatrix = TRUE)[[1]]
gse42133 <- getGEO("GSE42133", GSEMatrix = TRUE)[[1]]
gse111175 <- getGEO("GSE111175", GSEMatrix = TRUE)[[1]]
pheno_26415 <- pData(gse26415)
pheno_42133 <- pData(gse42133)
pheno_111175 <- pData(gse111175)
View(pheno_26415)
View(pheno_42133)
View(pheno_111175)

# GSE26415
labels_26415 <- data.frame(
  SampleID = rownames(pheno_26415),
  Label = trimws(pheno_26415$`sample type:ch1`)
)

# GSE42133
labels_42133 <- data.frame(
  SampleID = rownames(pheno_42133),
  Label = trimws(pheno_42133$`dx (diagnosis):ch1`)
)
labels_42133$Label <- tolower(labels_42133$Label)  # Convert "Control" to "control"

# GSE111175
labels_111175 <- data.frame(
  SampleID = rownames(pheno_111175),
  Label = trimws(pheno_111175$`diagnosis:ch1`)
)
labels_111175$Label <- ifelse(labels_111175$Label == "TD", "control", labels_111175$Label)
labels_26415$Batch <- "GSE26415"
labels_42133$Batch <- "GSE42133"
labels_111175$Batch <- "GSE111175"
all_labels <- bind_rows(labels_26415, labels_42133, labels_111175)
combined_exprs_labeled <- combined_exprs %>%
  left_join(all_labels, by = c("SampleID", "Batch"))
write.csv(combined_exprs_labeled, "C:\\Users\\GARIMA\\OneDrive\\Desktop\\dissertation\\ASD_Biomarker_discovery\\combined_expression_matrix_labeled.csv", row.names = FALSE)

# Convert all labels to lowercase
combined_exprs_labeled$Label <- tolower(combined_exprs_labeled$Label)

# Filter to keep only 'asd' and 'control' samples
combined_exprs_filtered <- combined_exprs_labeled %>%
  filter(Label %in% c("asd", "control"))

# Capitalize 'asd' to 'ASD' for readability
combined_exprs_filtered$Label <- ifelse(combined_exprs_filtered$Label == "asd", "ASD", "control")
write.csv(combined_exprs_filtered, "C:\\Users\\GARIMA\\OneDrive\\Desktop\\dissertation\\ASD_Biomarker_discovery\\combined_expression_matrix_final_label.csv", row.names = FALSE)
combined_exprs_filtered %>%
  count(Label)

#EDA
#1 Inspecting the data
# Dimensions
dim(combined_exprs_filtered)

# Preview the dataset
head(combined_exprs_filtered[, 1:6])  # first few genes and columns

# Check data types
str(combined_exprs_filtered)

# checking for missing values
# Count total missing values
sum(is.na(combined_exprs_filtered))
# Summary of expression values
summary(combined_exprs_filtered[, 1:5])  # first 5 genes only
#sample count per class
table(combined_exprs_filtered$Label)
#boxplot of expression distribution
set.seed(1)
genes_subset <- sample(colnames(combined_exprs_filtered)[1:(ncol(combined_exprs_filtered) - 3)], 100)
boxplot(combined_exprs_filtered[, genes_subset],
        las = 2,
        outline = FALSE,
        main = "Distribution of Expression for 100 Random Genes")

# Load required libraries
library(ggplot2)
library(dplyr)

# Read the labeled raw expression data
expr_data <- read.csv("C:/Users/GARIMA/OneDrive/Desktop/dissertation/ASD_Biomarker_discovery/combined_expression_matrix_final_label.csv")
# Extract expression matrix
expr_matrix <- expr_data %>%
  select(-SampleID, -Batch, -Label)
expr_matrix_log <- log2(expr_matrix + 1)
pca_result <- prcomp(expr_matrix_log, center = TRUE, scale. = TRUE)
# Combine PCA output with sample labels
pca_df <- as.data.frame(pca_result$x)
pca_df$Label <- expr_data$Label
ggplot(pca_df, aes(x = PC1, y = PC2, color = Label)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA of Gene Expression Data",
       x = "PC1",
       y = "PC2") +
  theme_minimal()
summary(pca_result)
pca_df$Batch <- expr_data$Batch  # Add Batch info

ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA Colored by Batch",
       x = "PC1", y = "PC2") +
  theme_minimal()

#batch correction
if (!requireNamespace("sva", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("sva")
}
library(sva)
library(dplyr)
# Read the final filtered labeled data
expr_data <- read.csv("C:/Users/GARIMA/OneDrive/Desktop/dissertation/ASD_Biomarker_discovery/combined_expression_matrix_final_label.csv")

# Select only expression columns
expr_matrix <- expr_data %>%
  select(-SampleID, -Batch, -Label)

# Log2-transform (if not already)
expr_matrix_log <- log2(expr_matrix + 1)

# Transpose so that rows = genes, columns = samples
expr_matrix_log_t <- t(expr_matrix_log)
batch <- expr_data$Batch
mod <- model.matrix(~ Label, data = expr_data) # Adjusting for condition
combat_expr <- ComBat(dat = expr_matrix_log_t, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)
combat_expr_df <- as.data.frame(t(combat_expr))
combat_expr_df$SampleID <- expr_data$SampleID
combat_expr_df$Batch <- expr_data$Batch
combat_expr_df$Label <- expr_data$Label

combat_expr_df <- as.data.frame(t(combat_expr))
combat_expr_df$SampleID <- expr_data$SampleID
combat_expr_df$Batch <- expr_data$Batch
combat_expr_df$Label <- expr_data$Label

# Save corrected data
write.csv(combat_expr_df, "batch_corrected_expression.csv", row.names = FALSE)

# PCA on corrected data
expr_matrix_corrected <- combat_expr_df %>% select(-SampleID, -Batch, -Label)
pca_corrected <- prcomp(expr_matrix_corrected, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca_corrected$x)
pca_df$Batch <- combat_expr_df$Batch
pca_df$Label <- combat_expr_df$Label

ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA After Batch Correction", x = "PC1", y = "PC2") +
  theme_minimal()

# post correction visualization
# Load libraries
library(ggplot2)

# PCA on corrected data
pca_corrected <- prcomp(combat_corrected_matrix, center = TRUE, scale. = TRUE)

# Prepare dataframe for plotting
pca_df_corrected <- as.data.frame(pca_corrected$x)
pca_df_corrected$SampleID <- rownames(pca_corrected$x)

# Merge with metadata
pca_df_corrected <- merge(pca_df_corrected, combined_exprs_filtered[, c("SampleID", "Label")], by = "SampleID")

# Plot
ggplot(pca_df_corrected, aes(x = PC1, y = PC2, color = Label)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA After Batch Correction (Colored by Label)",
       x = "PC1", y = "PC2") +
  theme_minimal()


# Load libraries
library(ggplot2)
library(dplyr)

# Load the corrected data
combat_expr_df <- read.csv("batch_corrected_expression.csv")

# Extract expression values only
expr_corrected <- combat_expr_df %>%
  select(-SampleID, -Batch, -Label)

# Run PCA
pca_result <- prcomp(expr_corrected, center = TRUE, scale. = TRUE)

# Create PCA dataframe
pca_df <- as.data.frame(pca_result$x)
pca_df$SampleID <- combat_expr_df$SampleID
pca_df$Label <- combat_expr_df$Label
pca_df$Batch <- combat_expr_df$Batch

# PCA Plot: Colored by Batch
ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA After Batch Correction",
       x = "PC1", y = "PC2") +
  theme_minimal()

# PCA Plot: Colored by Label (e.g., ASD vs Control)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Label)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA After Batch Correction (By Label)",
       x = "PC1", y = "PC2") +
  theme_minimal()

#DEG
install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR"))
# Load libraries
library(limma)
library(edgeR)

# Load your corrected expression file
combat_expr_df <- read.csv("batch_corrected_expression.csv")

# Separate metadata and expression matrix
expr_matrix <- combat_expr_df[, !(names(combat_expr_df) %in% c("SampleID", "Batch", "Label"))]
meta_data <- combat_expr_df[, c("SampleID", "Label")]
expr_matrix <- as.matrix(expr_matrix)
rownames(expr_matrix) <- meta_data$SampleID
# Convert Label to factor
meta_data$Label <- factor(meta_data$Label, levels = c("control", "ASD"))

# Design matrix for limma
design <- model.matrix(~ Label, data = meta_data)
# Fit linear model
fit <- lmFit(t(expr_matrix), design)
fit <- eBayes(fit)

# Get differentially expressed genes
results <- topTable(fit, coef = "LabelASD", number = Inf, adjust = "fdr")
# View top DE genes
head(results)

# Filter significant DE genes (adjusted p-value < 0.05 and |logFC| > 1)
de_genes <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
nrow(de_genes)
write.csv(results, "all_DE_results.csv", row.names = TRUE)
write.csv(de_genes, "significant_DE_genes.csv", row.names = TRUE)

de_genes_loose <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5, ]
nrow(de_genes_loose)
head(de_genes_loose[order(de_genes_loose$adj.P.Val), ])

library(ggplot2)

results$threshold <- as.factor(
  ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5, "Significant", "Not Significant")
)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 Adjusted P-value")
write.csv(de_genes_loose, "DEGs_filtered.csv", row.names = TRUE)

#functional enrichment
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2"))
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
gene_symbols <- rownames(de_genes_loose)  # Assuming rownames are gene symbols

entrez_ids <- bitr(gene_symbols,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene         = entrez_ids$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                ont          = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                readable     = TRUE)

# GO Dotplot
dotplot(ego, split = "ONTOLOGY") + facet_wrap(~ONTOLOGY)
