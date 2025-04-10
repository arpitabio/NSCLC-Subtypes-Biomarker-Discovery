if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks", force = TRUE)

# Create a dedicated directory
dir.create("TCGA_LUAD_Data")
setwd("TCGA_LUAD_Data")

# 1️⃣ Load required libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(Matrix)

# 2️⃣ Set batch size and total number of samples
batch_size <- 100  # Adjust as needed
total_samples <- 450
save_prefix <- "TCGA_LUAD_counts_batch"  # Prefix for saved files

# 3️⃣ Get full list of barcodes (sample IDs)
all_cases <- getResults(
  GDCquery(
    project = "TCGA-LUAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")  # Included sample type
  ), 
  cols = "cases"
)

# 4️⃣ Process data in batches
for (start in seq(1, total_samples, by = batch_size)) {
  
  # Define batch range
  end <- min(start + batch_size - 1, total_samples)
  batch_cases <- all_cases[start:end]
  
  # 5️⃣ Query for the current batch
  query <- GDCquery(
    project = "TCGA-LUAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal"),  # Included sample type
    barcode = batch_cases
  )
  
  # 6️⃣ Download in small chunks
  GDCdownload(query, files.per.chunk = 10)  # Prevents memory overload
  
  # 7️⃣ Prepare and save the batch data
  data <- GDCprepare(query)
  counts <- assay(data)  # Extract raw counts
  counts_sparse <- as(counts, "dgCMatrix")  # Convert to sparse matrix
  metadata <- colData(data)  # Extract clinical metadata
  
  # 8️⃣ Save batch to disk
  save_file <- paste0(save_prefix, "_", start, "_", end, ".rda")
  save(counts_sparse, metadata, file = save_file)
  
  # 9️⃣ Clear memory before next batch
  rm(data, counts, counts_sparse, metadata)
  gc()  # Garbage collection
}

# ✅ Done! Each batch is saved separately.

# 1️⃣ Load all batch files
batch_files <- list.files(pattern = "TCGA_LUAD_counts_batch_.*\\.rda$")
all_counts <- list()
all_metadata <- list()

for (file in batch_files) {
  load(file)  # Loads counts_sparse and metadata
  all_counts[[file]] <- counts_sparse
  all_metadata[[file]] <- metadata
}

# 2️⃣ Merge all batches into one dataset
final_counts <- do.call(cbind, all_counts) 
View(final_counts)
# Combine count matrices
final_metadata <- do.call(rbind, all_metadata)  # Combine metadata

# 3️⃣ Save final merged dataset
save(final_counts, final_metadata, file = "TCGA_LUAD_counts_FINAL.rda")

# ✅ Done! Use `load("TCGA_LUAD_counts_FINAL.rda")` to reload later.
# 1️⃣ Load the merged dataset
load("TCGA_LUAD_counts_FINAL.rda")  # Loads final_counts & final_metadata

# 2️⃣ Check available columns in metadata
colnames(final_metadata)

# 3️⃣ Count the number of Primary Tumor and Solid Tissue Normal samples
table(final_metadata$sample_type)

library(ggplot2)
ggplot(final_metadata, aes(x = sample_type)) + 
  geom_bar(fill = c("red", "blue")) + 
  ggtitle("Sample Distribution") + 
  xlab("Sample Type") + 
  ylab("Count")


rm(list = ls())  # Remove all variables
gc()  # Free up memory
load("TCGA_LUAD_counts_FINAL.rda")  # Load the dataset
ls()

library(Matrix)  # Load the Matrix package
counts <- as.matrix(final_counts)  # Convert to dense matrix (if needed)
dim(counts)  # Check dimensions
head(counts)  # View first few rows
head(final_metadata)  # View metadata
table(final_metadata$sample_type)  # Check counts of sample types

summary(final_counts)
any(is.na(final_counts))

rowSums(final_counts == 0)  # Genes with zero counts in all samples
summary(rowSums(final_counts))  # Summary of counts per gene
hist(log2(rowSums(final_counts) + 1), main="Log-transformed counts", xlab="Log2(counts+1)")


library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = final_counts, 
                              colData = final_metadata, 
                              design = ~ sample_type)
dds <- DESeq(dds)
res <- results(dds, contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"))
res
head(rownames(res))

library(stringr)

# Remove version suffix from Ensembl IDs
clean_ensembl_ids <- gsub("\\..*", "", rownames(res))

library(biomaRt)

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve gene symbols
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"), 
  filters = "ensembl_gene_id", 
  values = clean_ensembl_ids, 
  mart = mart
)


# Merge gene symbols with DE results
res$gene_symbol <- gene_annotations$hgnc_symbol[match(clean_ensembl_ids, gene_annotations$ensembl_gene_id)]

# Move gene_symbol to the first column
res <- res[, c("gene_symbol", colnames(res)[1:(ncol(res) - 1)])]

# View updated results
head(res)

# Remove rows where padj is NA
res_filtered <- res[!is.na(res$padj), ]

# Select upregulated genes
upregulated <- res_filtered[res_filtered$log2FoldChange > 1 & res_filtered$padj < 0.05, ]
rownames(upregulated) <- sub("\\..*", "", rownames(upregulated))

# Select downregulated genes
downregulated <- res_filtered[res_filtered$log2FoldChange < -1 & res_filtered$padj < 0.05, ]
rownames(downregulated) <- sub("\\..*", "", rownames(downregulated))

upregulated
downregulated

# Combine all significant DEGs
deg_list <- rbind(upregulated, downregulated)
dim(deg_list)
deg_list

# Convert rownames to a column for Upregulated genes
upregulated$Gene_Id <- rownames(upregulated)
rownames(upregulated) <- NULL  # Remove rownames
write.csv(upregulated, "Upregulated_genes_LUAD.csv", row.names = FALSE)

# Convert rownames to a column for Downregulated genes
downregulated$Gene_Id <- rownames(downregulated)
rownames(downregulated) <- NULL  # Remove rownames
write.csv(downregulated, "Downregulated_genes_LUAD.csv", row.names = FALSE)

# Save DEGs for future analysis
write.csv(deg_list, "TCGA_LUAD_DEGs.csv", row.names = FALSE)

# Quick summary
cat("Total DEGs:", nrow(deg_list), "\nUpregulated:", nrow(upregulated), "\nDownregulated:", nrow(downregulated))

#Volcano plot
library(ggplot2)

# Add labels for significant genes
res$Significance <- "Not Significant"
res$Significance[res$log2FoldChange > 1 & res$padj < 0.05] <- "Upregulated"
res$Significance[res$log2FoldChange < -1 & res$padj < 0.05] <- "Downregulated"

ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot of LUAD DEGs", x = "log2 Fold Change", y = "-log10 Adjusted p-value")

# Calculate the average expression (A) and log2 fold change (M)
A <- rowMeans(counts(dds))  # Average expression for each gene across samples
M <- res$log2FoldChange     # Log2 fold change for each gene

# MA plot: M vs A
plot(A, M, 
     pch = 20, 
     col = ifelse(res$padj < 0.05, "red", "gray"),  # Red for significant genes, gray for non-significant
     xlab = "Average Expression (A)", 
     ylab = "Log2 Fold Change (M)", 
     main = "MA Plot for LUAD DEGs (Primary Tumor vs Solid Tissue Normal)")
abline(h = 0, col = "black", lty = 2)  # Add a horizontal line at M = 0

# Load necessary libraries
library(DESeq2)
library(ggplot2)

# Perform PCA using transformed data (e.g., variance-stabilizing transformation)
vsd <- vst(dds, blind = TRUE)  # or rlog(dds, blind=TRUE)

# Get PCA data
pca_data <- plotPCA(vsd, intgroup = "sample_type", returnData = TRUE)

# Compute variance explained by each PC
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Plot PCA with variance explained
ggplot(pca_data, aes(x = PC1, y = PC2, color = sample_type)) +
  geom_point(size = 3) +
  labs(
    title = "PCA of TCGA-LUAD Data",
    x = paste0("PC1: ", percentVar[1], "%"),
    y = paste0("PC2: ", percentVar[2], "%")
  ) +
  theme_minimal()

library(clusterProfiler)
library(org.Hs.eg.db)

upregulated_genes <- rownames(res[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 0, ])
downregulated_genes <- rownames(res[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < 0, ])
# Remove version numbers
clean_upregulated <- sub("\\..*", "", upregulated_genes)
clean_downregulated <- sub("\\..*", "", downregulated_genes)

# Now, convert ENSEMBL to ENTREZ
up_entrez <- bitr(clean_upregulated, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_entrez <- bitr(clean_downregulated, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

head(clean_upregulated)
head(clean_downregulated)

failed_up <- clean_upregulated[!clean_upregulated %in% up_entrez$ENSEMBL]
failed_down <- clean_downregulated[!clean_downregulated %in% down_entrez$ENSEMBL]

print(failed_up)
print(failed_down)

library(biomaRt)

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Convert ENSEMBL to ENTREZ using biomaRt
up_entrez_biomart <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
                           filters = "ensembl_gene_id", 
                           values = failed_up, 
                           mart = mart)

down_entrez_biomart <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
                             filters = "ensembl_gene_id", 
                             values = failed_down, 
                             mart = mart)


# Combine successfully mapped ENTREZ IDs from both sources
final_up_entrez <- rbind(up_entrez, up_entrez_biomart)
final_down_entrez <- rbind(down_entrez, down_entrez_biomart)


colnames(up_entrez)  # Check column names from bitr
colnames(up_entrez_biomart)  # Check column names from biomaRt

colnames(down_entrez)
colnames(down_entrez_biomart)

colnames(up_entrez_biomart) <- c("ENSEMBL", "ENTREZID")
colnames(down_entrez_biomart) <- c("ENSEMBL", "ENTREZID")

# Remove duplicates (if any)
library(dplyr)
final_up_entrez <- rbind(up_entrez, up_entrez_biomart)
final_down_entrez <- rbind(down_entrez, down_entrez_biomart)

library(clusterProfiler)

# Gene Ontology (GO) enrichment for upregulated genes
go_up <- enrichGO(gene = final_up_entrez$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  keyType = "ENTREZID", 
                  ont = "BP", # Biological Process (You can also use "MF" for Molecular Function, "CC" for Cellular Component)
                  pAdjustMethod = "BH", 
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.05)

# Gene Ontology (GO) enrichment for downregulated genes
go_down <- enrichGO(gene = final_down_entrez$ENTREZID, 
                    OrgDb = org.Hs.eg.db, 
                    keyType = "ENTREZID", 
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.05)

# Visualize the top GO terms
dotplot(go_up, showCategory = 20) + ggtitle("GO Enrichment Biological Process - Upregulated Genes")
dotplot(go_down, showCategory = 20, label_format = 50) + ggtitle("GO Enrichment Biological Process - Downregulated Genes")

# Gene Ontology (GO) enrichment for upregulated genes
go_up_MF <- enrichGO(gene = final_up_entrez$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  keyType = "ENTREZID", 
                  ont = "MF", # Biological Process (You can also use "MF" for Molecular Function, "CC" for Cellular Component)
                  pAdjustMethod = "BH", 
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.05)

# Gene Ontology (GO) enrichment for downregulated genes
go_down_MF <- enrichGO(gene = final_down_entrez$ENTREZID, 
                    OrgDb = org.Hs.eg.db, 
                    keyType = "ENTREZID", 
                    ont = "MF", 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.05)

# Visualize the top GO terms
dotplot(go_up_MF, showCategory = 20, label_format = 50) + ggtitle("GO Enrichment Molecular Function - Upregulated Genes")
dotplot(go_down_MF, showCategory = 20, label_format = 50) + ggtitle("GO Enrichment Molecular Function - Downregulated Genes")



# Gene Ontology (GO) enrichment for upregulated genes
go_up_CC <- enrichGO(gene = final_up_entrez$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  keyType = "ENTREZID", 
                  ont = "CC", # Biological Process (You can also use "MF" for Molecular Function, "CC" for Cellular Component)
                  pAdjustMethod = "BH", 
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.05)

# Gene Ontology (GO) enrichment for downregulated genes
go_down_CC <- enrichGO(gene = final_down_entrez$ENTREZID, 
                    OrgDb = org.Hs.eg.db, 
                    keyType = "ENTREZID", 
                    ont = "CC", 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.05)


# Visualize the top GO terms
dotplot(go_up_CC, showCategory = 20, label_format = 50) + ggtitle("GO Enrichment Cellular Component - Upregulated Genes")
dotplot(go_down_CC, showCategory = 20, label_format = 50) + ggtitle("GO Enrichment Cellular Component - Downregulated Genes")



# KEGG pathway enrichment for upregulated genes
kegg_up <- enrichKEGG(gene = final_up_entrez$ENTREZID, 
                      organism = "hsa", 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 0.05)

# KEGG pathway enrichment for downregulated genes
kegg_down <- enrichKEGG(gene = final_down_entrez$ENTREZID, 
                        organism = "hsa", 
                        pAdjustMethod = "BH", 
                        pvalueCutoff = 0.05)

# Visualize KEGG pathways
# KEGG Pathway Enrichment - Upregulated Genes (Green-Blue)
barplot(kegg_up, showCategory = 15, label_format = 50) + 
  ggtitle("KEGG Pathway Enrichment - Upregulated Genes") +
  scale_fill_gradient(low = "green", high = "blue")  # Green to Blue gradient

# KEGG Pathway Enrichment - Downregulated Genes (Yellow-Violet)
barplot(kegg_down, showCategory = 15, label_format = 50) + 
  ggtitle("KEGG Pathway Enrichment - Downregulated Genes") +
  scale_fill_gradient(low = "orange", high = "violet")  # Yellow to Violet gradient


# Ensure deg_list is in data frame format
# Convert DESeqResults to a data frame
deg_df <- as.data.frame(deg_list)

# Ensure ENSEMBL IDs are in a column (currently in row names)
deg_df$ENSEMBL <- rownames(deg_df)

# Combine the ENTREZ mappings from upregulated and downregulated genes
entrez_mapping <- rbind(final_up_entrez, final_down_entrez)

# Merge with ENTREZ mapping (keep only mapped genes)
merged_degs <- merge(deg_df, entrez_mapping, by = "ENSEMBL")

# Select only ENTREZ ID and log2FoldChange
ranked_genes <- merged_degs[, c("ENTREZID", "log2FoldChange")]

# Remove NAs (if any ENTREZ IDs failed to map)
ranked_genes <- na.omit(ranked_genes)

# Convert to named numeric vector
ranked_genes_vector <- setNames(ranked_genes$log2FoldChange, ranked_genes$ENTREZID)

# Sort in decreasing order (GSEA requires ranked list)
ranked_genes_vector <- sort(ranked_genes_vector, decreasing = TRUE)

library(clusterProfiler)
library(org.Hs.eg.db)


library(dplyr)

# Convert to data frame, group by ENTREZ ID, and take the mean log2FC
ranked_genes_df <- data.frame(ENTREZID = names(ranked_genes_vector), logFC = ranked_genes_vector)

ranked_genes_unique <- ranked_genes_df %>%
  group_by(ENTREZID) %>%
  summarize(logFC = mean(logFC, na.rm = TRUE)) %>%
  ungroup()

# Convert back to named vector
ranked_genes_vector <- setNames(ranked_genes_unique$logFC, ranked_genes_unique$ENTREZID)

sum(duplicated(names(ranked_genes_vector)))

ranked_genes_vector <- sort(ranked_genes_vector, decreasing = TRUE)

head(ranked_genes_vector)  # Should show highest log2FC first
tail(ranked_genes_vector)  # Should show lowest log2FC last

str(ranked_genes_vector)

#GSEA Pathway Enrichment

gsea_go <- gseGO(
  geneList = ranked_genes_vector,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)


# Visualize results
dotplot(gsea_go, showCategory = 20, label_format = 50) + 
  ggtitle("GSEA - GO Biological Process")

dim(as.data.frame(gsea_go))

library(ggplot2)

# Convert GSEA results to a dataframe
gsea_df <- as.data.frame(gsea_go)

# Check if gsea_df has data
head(gsea_df)

ggplot(gsea_df[1:20, ], aes(x = reorder(Description, NES), y = NES, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") + 
  labs(title = "GSEA - 20 GO Biological Process", x = "GO Terms", y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),           # Increase overall text size
    axis.text.x = element_text(size = 12),    # Increase x-axis text size
    axis.text.y = element_text(size = 12),    # Increase y-axis text size
    axis.title = element_text(size = 14, face = "bold"),  # Bold axis labels
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5) # Centered bold title
  )


library(enrichplot)

gsea_go_sim <- pairwise_termsim(gsea_go)  # Calculate similarity
emapplot(gsea_go_sim, showCategory = 20)  # Now plot

gseaplot2(gsea_go, geneSetID = 1, title = gsea_go$Description[1])

gsea_kegg <- gseKEGG(geneList = ranked_genes_vector, 
                     organism = "hsa", 
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     pvalueCutoff = 0.05, 
                     pAdjustMethod = "BH")

dotplot(gsea_kegg, showCategory = 20) + ggtitle("GSEA - KEGG Pathways")


head(as.data.frame(gsea_go))
write.csv(as.data.frame(gsea_go), "GSEA_GO_Results.csv", row.names = FALSE)
write.csv(as.data.frame(gsea_kegg), "GSEA_KEGG_Results.csv", row.names = FALSE)

dim(as.data.frame(gsea_kegg))


library(ggplot2)

# Convert GSEA results to a data frame
gsea_kegg_df <- as.data.frame(gsea_kegg)

# Select the top 20 enriched pathways
top_kegg <- gsea_kegg_df[1:20, ]

# Create the bar plot
ggplot(top_kegg, aes(x = reorder(Description, NES), y = NES, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip for better readability
  scale_fill_gradient(low = "blue", high = "red") + 
  labs(title = "GSEA - 20 KEGG Pathway Enrichment", x = "KEGG Pathways", y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),          # Increase overall text size
    axis.text.y = element_text(size = 12),   # Increase y-axis label size
    axis.text.x = element_text(size = 12),   # Increase x-axis label size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 12)    # Increase legend text size
  )

library(enrichplot)
emapplot(pairwise_termsim(gsea_kegg), showCategory = 20)

ridgeplot(gsea_kegg, showCategory = 20) + 
  ggtitle("Distribution of Genes in Enriched KEGG Pathways") +
  theme(
    axis.text.y = element_text(size = 9.5),  # Increase y-axis text size
    axis.text.x = element_text(size = 9.5),  # Increase x-axis text size
    plot.title = element_text(size = 14, face = "bold")  # Increase title size
  )


# Load required package
library(org.Hs.eg.db)

# Get all valid Entrez Gene IDs
valid_entrez_list <- keys(org.Hs.eg.db, keytype = "ENTREZID")

# Convert to numeric and remove long IDs (>6 digits)
valid_entrez_genes <- names(ranked_genes_vector)
valid_entrez_genes <- valid_entrez_genes[nchar(valid_entrez_genes) <= 6]  # Remove long ones

# Keep only genes that exist in Entrez database
valid_entrez_genes <- valid_entrez_genes[valid_entrez_genes %in% valid_entrez_list]

# Subset the ranked vector
ranked_genes_vector_cleaned <- ranked_genes_vector[valid_entrez_genes]

# Check if filtering worked
head(ranked_genes_vector_cleaned)


length(unique(names(ranked_genes_vector_cleaned)))  # Should be > 0

if (!require("DOSE")) BiocManager::install("DOSE")
library(DOSE)

gsea_do <- enrichDO(gene = names(ranked_genes_vector_cleaned), pvalueCutoff = 0.05)
head(gsea_do@result)

gsea_kegg0 <- enrichKEGG(gene = names(ranked_genes_vector_cleaned), organism = "hsa")
head(gsea_kegg0@result)

dim(gsea_do@result)  # Check the number of enriched disease terms
dim(gsea_kegg0@result)  # Check the number of enriched pathways

head(gsea_do@result[, c("ID", "Description", "p.adjust")])  # Top diseases
head(gsea_kegg0@result[, c("ID", "Description", "p.adjust")])  # Top pathways

#Barplot of Top 10 Enriched Diseases
barplot(gsea_do, showCategory = 10)

#Barplot of Top 10 KEGG Pathways
barplot(gsea_kegg0, showCategory = 10)

lung_cancer_do <- subset(gsea_do@result, grepl("lung cancer|lung carcinoma", Description, ignore.case = TRUE))
print(lung_cancer_do)

lung_cancer_kegg <- subset(gsea_kegg0@result, grepl("lung cancer|NSCLC|SCLC|carcinoma", Description, ignore.case = TRUE))
print(lung_cancer_kegg)

lung_cancer_genes_do <- unique(unlist(strsplit(lung_cancer_do$geneID, "/")))
print(lung_cancer_genes_do)

lung_cancer_genes_kegg <- unique(unlist(strsplit(lung_cancer_kegg$geneID, "/")))
print(lung_cancer_genes_kegg)

ego <- enrichGO(gene = lung_cancer_genes_do, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH")
dotplot(ego, showCategory = 10)

#Convert DO-Enriched Genes
lung_cancer_gene_symbols_do <- mapIds(org.Hs.eg.db, 
                                      keys = lung_cancer_genes_do, 
                                      column = "SYMBOL", 
                                      keytype = "ENTREZID", 
                                      multiVals = "first")
print(lung_cancer_gene_symbols_do)
# Convert lung_cancer_gene_symbols_do into a dataframe
lung_cancer_df <- data.frame(
  Entrez_ID = names(lung_cancer_gene_symbols_do),  # Extract Entrez IDs
  Gene_Symbol = as.character(lung_cancer_gene_symbols_do),  # Extract Gene Symbols
  stringsAsFactors = FALSE
)

# Display the dataframe
print(lung_cancer_df)

ego_2 <- enrichGO(gene = lung_cancer_genes_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH")
dotplot(ego_2, showCategory = 10)

#Convert DO-Enriched Genes
lung_cancer_gene_symbols_kegg <- mapIds(org.Hs.eg.db, 
                                      keys = lung_cancer_genes_kegg, 
                                      column = "SYMBOL", 
                                      keytype = "ENTREZID", 
                                      multiVals = "first")
print(lung_cancer_gene_symbols_kegg)
# Convert lung_cancer_gene_symbols_do into a dataframe
lung_cancer_kegg_df <- data.frame(
  Entrez_ID = names(lung_cancer_gene_symbols_kegg),  # Extract Entrez IDs
  Gene_Symbol = as.character(lung_cancer_gene_symbols_kegg),  # Extract Gene Symbols
  stringsAsFactors = FALSE
)

# Display the dataframe
print(lung_cancer_kegg_df)



kegg_enrichment <- enrichKEGG(gene = lung_cancer_genes_kegg, 
                              organism = "hsa", 
                              pAdjustMethod = "BH") 

dotplot(kegg_enrichment, showCategory = 10)

write.csv(lung_cancer_kegg_df, "lung_cancer_kegg_biomarkers.csv", row.names = FALSE)
write.csv(lung_cancer_df, "lung_cancer_biomarkers.csv", row.names = FALSE)

# Find common Entrez IDs
common_entrez_ids <- intersect(lung_cancer_kegg_df[,1], lung_cancer_df[,1])

# Extract rows where Entrez ID is in the common list
final_biomarkers_df <- lung_cancer_kegg_df[lung_cancer_kegg_df[,1] %in% common_entrez_ids, ]

# Rename columns for clarity (assuming 1st column is Entrez ID and 2nd is Gene Symbol)
colnames(final_biomarkers_df) <- c("Entrez_ID", "Gene_Symbol")

# Save the final biomarkers to CSV
write.csv(final_biomarkers_df, "final_biomarkers.csv", row.names = FALSE)

# Print the final common biomarkers
print(final_biomarkers_df)

