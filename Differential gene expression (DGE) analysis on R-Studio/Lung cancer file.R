library(DESeq2)
library(ggplot2)
library(pheatmap)

# Read the count matrix
count_data <- read.csv("Lung_cancer_dataset.csv", row.names = 1)
head(count_data)

# Sample metadata
sample_info <- data.frame(
  row.names = colnames(count_data),
  condition = c("Knockdown", "Knockdown", "Knockdown", "WildType", "WildType", "WildType")
)

#Create DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~ condition)

dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- results(dds)


# Order by adjusted p-value and save
res_ordered <- res[order(res$padj), ]
write.csv(as.data.frame(res_ordered), file = "DEG_results.csv")


#PCA Plot
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

#MA Plot
plotMA(res, ylim = c(-5, 5))

#Volcano Plot
res_df <- as.data.frame(res)
res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant)) +
  scale_color_manual(values = c("black", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "log2(Fold Change)", y = "-log10(padj)")


#Heatmap of Top 30 DEGs
topgenes <- head(order(res$padj), 30)
mat <- assay(vsd)[topgenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = sample_info)

upregulated <- subset(res_ordered, padj < 0.05 & log2FoldChange > 1)
write.csv(as.data.frame(upregulated), file = "Upregulated_Genes.csv")

downregulated <- subset(res_ordered, padj < 0.05 & log2FoldChange < -1)
write.csv(as.data.frame(downregulated), file = "Downregulated_Genes.csv")







#GO and KEGG Enrichment
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

#Prepare Gene List
# Load upregulated genes (example)
up_genes <- read.csv("upregulated_genes.csv", row.names = 1)

# Extract Ensembl IDs and convert to Entrez IDs
ensembl_ids <- rownames(up_genes)
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                       column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")

# Remove NAs
gene_symbols <- na.omit(gene_symbols)

#GO Enrichment (Biological Process)
ego <- enrichGO(gene         = gene_symbols,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# View results
head(ego)

# Plot top GO terms
barplot(ego, showCategory = 15, title = "GO Enrichment - Biological Process")

#KEGG Pathway Enrichment
ekegg <- enrichKEGG(gene         = gene_symbols,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05)

# View results
head(ekegg)

# KEGG dotplot
dotplot(ekegg, showCategory = 15, title = "KEGG Pathway Enrichment")

write.csv(as.data.frame(ego), "GO_Enrichment.csv")
write.csv(as.data.frame(ekegg), "KEGG_Enrichment.csv")




# Total DEGs
deg_all <- subset(res_ordered , padj < 0.05)

# Downregulated
deg_down <- subset(deg_all, log2FoldChange < -1)

# Upregulated
deg_up <- subset(deg_all, log2FoldChange > 1)


# Number of DEGs
cat("Total significant DEGs (padj < 0.05):", nrow(deg_all), "\n")

# Downregulated genes
cat("Downregulated genes (log2FC < -1):", nrow(deg_down), "\n")

# Upregulated genes
cat("Upregulated genes (log2FC > 1):", nrow(deg_up), "\n")
head(deg_down)
head(deg_up)




# Connect to Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")  # For human

# Extract unique Ensembl IDs from your results
gene_ids <- rownames(res_ordered)

# Get mapping: Ensembl ID → Gene Symbol
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = ensembl
)

# Add gene IDs as a column in your DEG tables
deg_down$ensembl_gene_id <- rownames(deg_down)
deg_up$ensembl_gene_id <- rownames(deg_up)


# Convert to data frame before merging
deg_down_df <- as.data.frame(deg_down)
deg_down_df$ensembl_gene_id <- rownames(deg_down_df)

deg_up_df <- as.data.frame(deg_up)
deg_up_df$ensembl_gene_id <- rownames(deg_up_df)

# Now merge with gene_info
deg_down_annot <- merge(deg_down_df, gene_info, by = "ensembl_gene_id", all.x = TRUE)
deg_up_annot <- merge(deg_up_df, gene_info, by = "ensembl_gene_id", all.x = TRUE)



write.csv(deg_down_annot, "downregulated_genes_with_names.csv", row.names = FALSE)
write.csv(deg_up_annot, "upregulated_genes_with_names.csv", row.names = FALSE)
