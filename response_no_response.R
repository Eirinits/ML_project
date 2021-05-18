library(tidyverse)
require(tibble)
require(ggplot2)
library(umap)
# library(stringr)
library(VennDiagram)
library(RColorBrewer)
# library("preprocessCore")
# library(VIM)
# library(FactoMineR)
library(DESeq2)
library(pheatmap)
library(PCAtools)


metadata <- metadata_GSE106992 %>% 
  dplyr::select(refinebio_accession_code,time,treatment,description,pasi75resp) %>% 
  mutate(responder = factor(paste0(description,"_",treatment,"_",pasi75resp))) 

metadata$responder <- sub(" Week ","", metadata$responder)
metadata$responder <- sub("Etanercept","Etan", metadata$responder)
metadata$responder <- sub("Ustekinumab 90 mg","Ust", metadata$responder)

# psoriasis_samples <- c(metadata$refinebio_accession_code,"Gene")
# GSE106992 <- subset(GSE106992, select = psoriasis_samples)

# Make the data in the order of the metadata
#expression_df <- data_filtered %>%
#  column_to_rownames("Gene") %>% 
#  dplyr::select(metadata$refinebio_accession_code)
expression_df <- data_filtered
# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)  

gene_matrix <- round(filtered_expression_df)
gene_matrix <- gene_matrix[apply(gene_matrix[,-1], 1, function(x) !all(x==0)),]
gene_matrix[gene_matrix <0] <- 0

ddset <- DESeqDataSetFromMatrix(
  countData = gene_matrix,
  colData = metadata,
  design = ~ responder 
)

dds_norm <- vst(ddset)

normalized_counts <- assay(dds_norm) %>%
  t()

umap_results <- umap::umap(normalized_counts)

# Make into data frame for plotting with `ggplot2`
# The UMAP values we need for plotting are stored in the `layout` element
umap_plot_df <- data.frame(umap_results$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("refinebio_accession_code") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(metadata, by = "refinebio_accession_code")

umap_plot_df
# Plot using `ggplot()` function
ggplot(
  umap_plot_df,
  aes(
    x = X1,
    y = X2,
    color = responder# label points with different colors for each `subgroup`
  )
) +
  geom_point() # This tells R that we want a scatterplot Plot individual points to make a scatterplot

dds <- DESeq(ddset)

unique(metadata$responder)
deseq_results_UST_Y <- results(dds, 
                               contrast = c("responder",
                                            "LS0_Ust_Yes",
                                            "NL0_Ust_Yes"))

deseq_df_UST_Y <- deseq_results_UST_Y %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::mutate(threshold = padj < 0.01) %>%
  dplyr::arrange(dplyr::desc(log2FoldChange)) %>%
  filter(threshold == "TRUE") %>%  
  filter(abs(log2FoldChange) > 1)

volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df_UST_Y,
  lab = deseq_df_UST_Y$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)

volcano_plot

deseq_results_UST_N <- results(dds, contrast = c("responder",
                                                 "LS0_Ust_No",
                                                 "NL0_Ust_No"))

deseq_df_UST_N <- deseq_results_UST_N %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::mutate(threshold = padj < 0.01) %>%
  dplyr::arrange(dplyr::desc(log2FoldChange)) %>%
  filter(threshold == "TRUE") %>%  
  filter(abs(log2FoldChange) > 1)

volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df_UST_N,
  lab = deseq_df_UST_N$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)

volcano_plot

deseq_results_Et_N <- results(dds, contrast = c("responder",
                                                 "LS0_Etan_No",
                                                 "NL0_Etan_No"))

deseq_df_Et_N <- deseq_results_Et_N %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::mutate(threshold = padj < 0.01) %>%
  dplyr::arrange(dplyr::desc(log2FoldChange)) %>%
  filter(threshold == "TRUE") %>%  
  filter(abs(log2FoldChange) > 1)

volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df_Et_N,
  lab = deseq_df_Et_N$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)

volcano_plot

deseq_results_Et_Y <- results(dds, contrast = c("responder",
                                                "LS0_Etan_Yes",
                                                "NL0_Etan_Yes"))

deseq_df_Et_Y <- deseq_results_Et_Y %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::mutate(threshold = padj < 0.01) %>%
  dplyr::arrange(dplyr::desc(log2FoldChange)) %>%
  filter(threshold == "TRUE") %>%  
  filter(abs(log2FoldChange) > 1)

volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df_Et_Y,
  lab = deseq_df_Et_Y$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)

volcano_plot

deseq_results_NonResp <- results(dds, contrast = c("responder",
                                                "LS0_Etan_Yes",
                                                "LS0_Etan_No"))

deseq_df_NonResp <- deseq_results_NonResp %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::mutate(threshold = padj < 0.05) %>%
  dplyr::arrange(dplyr::desc(log2FoldChange)) %>%
  filter(threshold == "TRUE") 

volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df_Et_Y,
  lab = deseq_df_Et_Y$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)

volcano_plot

library(gplots)
venn(
  data = list(deseq_df_Et_Y$Gene,deseq_df_Et_N$Gene,deseq_df_UST_Y$Gene,deseq_df_UST_N$Gene),
  #category.names = c("Etan_responders","Etan_Nonresponders","Ust_responders","Ust_Nonresponders"),
# filename = "DEGs.png",
#  output = TRUE
)

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

diff_ust <- outersect(deseq_df_UST_Y$Gene,deseq_df_UST_N$Gene)
diff_et <- outersect(deseq_df_Et_Y$Gene,deseq_df_Et_N$Gene)
comm_no <- intersect(deseq_df_Et_N$Gene,deseq_df_UST_N$Gene)
comm_yes <- intersect(deseq_df_Et_Y$Gene,deseq_df_UST_Y$Gene)
setdiff(comm_yes,comm_no)

# Pca ---------------------------------------------------------------------

target <- colnames(data_filtered)
metadata_ordered <- metadata[match(target, metadata$refinebio_accession_code),]
metadata_ordered <- metadata %>% 
  column_to_rownames("refinebio_accession_code")

data <- data_filtered[,rownames(metadata_ordered)]
p <- pca(data_filtered, metadata = metadata_ordered, removeVar = 0.2, scale = FALSE)

screeplot(p, axisLabSize = 18, titleLabSize = 22)
pairsplot(p,colby = "pasi75resp")

biplot(p,
       legendPosition = 'bottom',
       colby = 'responder',
       lab = NULL,
       shape = 'pasi75resp')

pairsplot(p, components = getComponents(p, c(1:3)),
          triangle = TRUE, trianglelabSize = 12,
          hline = 0, vline = 0, pointSize = 2,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'responder', shape = "pasi75resp",  plotaxes = FALSE,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'),
          legendPosition = "bottom")

plotloadings(p, rangeRetain = 0.01,
             labSize = 4.0, title = 'Loadings plot',
             subtitle = 'PC1, PC2, PC3, PC4, PC5', caption = 'Top 1% variables',
             shape = 24, col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)