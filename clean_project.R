##################################### FIGURE 5 #################################
library(limma)
library(ggplot2)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(biomaRt)
library(edgeR) 

setwd("~/git/R")


metadata <- read.table("datasets/updated_metadata.txt", header = TRUE, sep = "")
counts <- read.table("datasets/updated_counts.txt", header = TRUE, sep = " ", row.names = 1)

# Add a new column in metadata to identify responders and non-responders
metadata$response <- ifelse(metadata$pdx %in% c("T111", "T113", "PL-015"), "responder", "non_responder")


metadata$group <- interaction(metadata$response, metadata$treatment)
design <- model.matrix(~ 0 + metadata$group)
colnames(design) <- make.names(colnames(design))

# Create a DGEList 
dge <- DGEList(counts = counts)

# Normalize the DGEList object using voom
v <- voom(dge, design, plot=TRUE)

# Fit the linear model using the voom normalized data
fit <- lmFit(v, design)

# Update the contrast matrix based on the column names of the design matrix
contrast_matrix <- makeContrasts(P4_responders = metadata.groupresponder.P4 - metadata.groupresponder.CTRL,
                                 P4_nonresponders = metadata.groupnon_responder.P4 - metadata.groupnon_responder.CTRL,
                                 levels = design)

# Apply the contrasts and perform eBayes
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract results for responders and nonresponders
responders <- topTable(fit2, coef = "P4_responders", n = Inf)
nonresponders <- topTable(fit2, coef = "P4_nonresponders", n = Inf)

# Filter topTable results based on log2(fold change) and adjusted p-value
filtered_responders <- responders[responders$logFC > -10 & responders$logFC < 10,]
filtered_nonresponders <- nonresponders[nonresponders$logFC > -10 & nonresponders$logFC < 10, ]

# Convert ENSEMBL gene IDs to gene symbols
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
responder_genesymbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id", values = rownames(filtered_responders), mart = ensembl)
nonresponder_genesymbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id", values = rownames(filtered_nonresponders), mart = ensembl)

# Merge gene symbols with filtered results
filtered_responders <- merge(responder_genesymbols, filtered_responders, by.x = "ensembl_gene_id", by.y = "row.names")
filtered_nonresponders <- merge(nonresponder_genesymbols, filtered_nonresponders, by.x = "ensembl_gene_id", by.y = "row.names")


# Create a table of the duplicated names
dupes <- table(filtered_responders$external_gene_name)
dupes <- names(dupes[dupes > 1])

# Add a suffix to the duplicate names
for (dupe in dupes) {
  indices <- which(filtered_responders$external_gene_name == dupe)
  filtered_responders$external_gene_name[indices] <- paste0(dupe, "_", seq_along(indices))
}

rownames(filtered_responders) <- filtered_responders$external_gene_name

# Same for non responders
dupes <- table(filtered_nonresponders$external_gene_name)
dupes <- names(dupes[dupes > 1])
for (dupe in dupes) {
  indices <- which(filtered_nonresponders$external_gene_name == dupe)
  filtered_nonresponders$external_gene_name[indices] <- paste0(dupe, "_", seq_along(indices))
}
rownames(filtered_nonresponders) <- filtered_nonresponders$external_gene_name


# Define the color mapping based on logFC and P-value
keyvals.colour <- ifelse(filtered_responders$logFC > 0.5 & filtered_responders$P.Value < 0.05, 'red',
                         ifelse(filtered_responders$logFC < -0.5 & filtered_responders$P.Value < 0.05, 'blue',
                                'gray'))


# Replace NA values with gray
keyvals.colour[is.na(keyvals.colour)] <- 'gray'
names(keyvals.colour)[keyvals.colour == 'red'] <- 'log2(FC) > 0.5'  
#names(keyvals.colour)[keyvals.colour == 'gray'] <- '-'  
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'log2(FC) < 0.5' 

# a
volcano_responders <- EnhancedVolcano(filtered_responders,
                                      lab = rownames(filtered_responders),
                                      x = 'logFC',
                                      y = 'P.Value',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      title = "P4-induced changes in responders",
                                      xlab = "Log2(fold change)",
                                      ylab = "-log10(P-value)",
                                      colCustom = keyvals.colour)
print(volcano_responders)

# b
volcano_nonresponders <- EnhancedVolcano(filtered_nonresponders,
                                         lab = rownames(filtered_nonresponders),
                                         x = 'logFC',
                                         y = 'P.Value',
                                         pCutoff = 0.05,
                                         FCcutoff = 0.5,
                                         title = "P4-induced changes in non-responders",
                                         xlab = "Log2(fold change)",
                                         ylab = "-log10(P-value)",
                                         colCustom = keyvals.colour)

print(volcano_nonresponders)

# c
msigdb <- msigdbr(species = "Homo sapiens", category = "H")
gene_list <- filtered_responders$logFC
names(gene_list) <- rownames(filtered_responders) 

# Data frame for the term-to-gene mapping
term2gene <- data.frame(term = msigdb$gs_name, gene = msigdb$gene_symbol)

# Data frame for the term-to-name mapping
term2name <- data.frame(term = msigdb$gs_name, name = msigdb$gs_name)

# Perform gene set enrichment analysis using GSEA
sorted_gene_list <- sort(gene_list, decreasing = TRUE)
gsea_responders <- GSEA(geneList = sorted_gene_list, TERM2GENE = term2gene, TERM2NAME = term2name)

# Get the enriched pathways
enriched_pathways <- gsea_responders$ID[gsea_responders$p.adjust < 0.05]

# Create a dot plot for the enriched pathways
enriched_df <- gsea_responders[order(gsea_responders$NES), ]
size <- ifelse(enriched_df$setSize < 50, 1,
               ifelse(enriched_df$setSize < 100, 2, 3))

ggplot(enriched_df, aes(x = NES, y = ID, color = p.adjust, size = setSize)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point() +
  scale_color_continuous(name = "p.adjust", breaks = c(0.01, 0.02, 0.03, 0.04), labels = c("0.01", "0.02", "0.03", "0.04")) +
  scale_size_continuous(name = "SetSize", breaks = c(50, 100, 150), labels = c("50", "100", "150")) +
  labs(title = "Enriched pathways in P4 vs CTRL (in Responders)", x = "Normalized Enrichment Score", y = "Pathway") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8, hjust = 0),
        axis.text.x = element_text(size = 8),
        legend.position = "right",
        legend.text = element_text(size = 8)) +
  guides(color = guide_colorbar(title = "p.adjust", label = TRUE),
         size = guide_legend(title = "SetSize"))

