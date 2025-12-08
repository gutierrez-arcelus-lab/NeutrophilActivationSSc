#!/usr/bin/env Rscript

###############################################################################
# SYSTEMIC SCLEROSIS — GSEA + DEG VISUALIZATION PIPELINE
# Author: Gonzalo Villanueva-Martin
# Date: 2025
# Description:
#   This script aggregates results from multiple SSc datasets, performs:
#     - Neutrophil/platelet GSEA extraction
#     - NES visualization across datasets
#     - Tissue annotation
#     - DEG extraction for neutrophil/platelet genes
#     - Annotated heatmaps with p-value symbols
#     - Enrichment plots panel (2x4)
###############################################################################


################################################################################
# Load required packages
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(reshape2)
  library(tibble)
  library(openxlsx)
  library(GhibliBrewer)
  library(dplyr)
  library(purrr)
  library(pheatmap)
  library(readr)
  library(stringr)
  library(patchwork)
  library(fgsea)
  library(msigdbr)
  library(cowplot)
})


################################################################################
# FIGURE COLORS + DATASET ORDERING
################################################################################

# Colorblind-friendly palette for tissues
tissue_cols <- c(
  "Skin"  = "#5E3C99",   # Purple – Skin
  "Blood" = "#D73027"   # Red – Blood
)

# Desired order for dataset appearance across all figures
tissue_order <- c(
  "GSE130953","GSE231691", # Blood
  "GSE249550", "GSE181549", "GSE130955","GSE95065" # Skin
)



################################################################################
# LOAD GSEA RESULTS
################################################################################

gsea_res <- list.files("results/fgsea", pattern = "\\.txt$", full.names = TRUE) %>%
  setNames(gsub("REACTOME_results_|\\.txt$", "", basename(.))) %>%
  lapply(read_tsv)

gsea_res <- gsea_res[intersect(names(gsea_res), tissue_order)]

################################################################################
# LOAD NEUTROPHIL GENESET (MSigDB REACTOME)
################################################################################

neutroReacPath <- read_tsv("REACTOME_NEUTROPHIL_DEGRANULATION.v2025.1.Hs.tsv")
genesNeu <- str_split(neutroReacPath[17, 2], ",")[[1]]


################################################################################
# EXTRACT NEUTROPHIL/PLATELET RELATED GSEA RESULTS
################################################################################

neutro_list <- imap(gsea_res, ~{
  df <- .x
  name <- .y
  
  neu_term <- grep("(NEUTROPHIL|PLATELET)", df$pathway,
                   ignore.case = TRUE, value = TRUE)
  
  df %>%
    filter(pathway %in% neu_term) %>%
    filter(padj < 0.1) %>%
    arrange(desc(NES)) %>%
    mutate(
      dataset = name,
      Tissue = case_when(
        grepl("GSE95065|GSE181549|GSE130955|GSE249550", dataset) ~ "Skin",
        grepl("GSE130953|GSE231691", dataset) ~ "Blood",
        TRUE ~ "Lung"
      )
    )
})

# Combine all datasets
neutro_all <- bind_rows(neutro_list)
plate_df <- neutro_all %>% filter(pathway != "REACTOME_NEUTROPHIL_DEGRANULATION")


################################################################################
# EXTRACT ONLY NEUTROPHIL DEGRANULATION PATHWAY
################################################################################

neu_term <- grep("NEUTROPHIL", neutro_all$pathway,
                 ignore.case = TRUE, value = TRUE)

neu_only <- neutro_all %>%
  filter(pathway == neu_term) %>%
  arrange(desc(padj))

# Compute ratio of detected leading-edge genes
neu_only$ratio <- neu_only$size / length(genesNeu)


################################################################################
# SIDE ANNOTATION DF
################################################################################

anno_df <- neu_only %>%
  select(dataset, Tissue, NES)

anno_pla_df <- plate_df %>%
  select(dataset, Tissue, NES)

###############################################################################
# MAIN NES BARPLOT (NEUTROPHIL DEGRANULATION)
###############################################################################

NESplot_2 <- ggplot(neu_only, aes(factor(dataset, levels = c(rev(tissue_order))), NES)) +
  geom_col(aes(fill = Tissue),
           color = "black", alpha = 0.8) +
  scale_fill_manual(
    values = tissue_cols,
    name = "NES FDR 10%"
  ) +
  coord_flip() +
  labs(
    x = "",
    y = "Normalized Enrichment Score",
    title = "Reactome Neutrophil Degranulation\nacross SSc datasets"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

ggsave("results/plots4publication/NeuDegranGSEAres_v2.png",
       plot = NESplot_2, dpi = 300, width = 6, height = 5.5)

###############################################################################
# DOTPLOT — ALL NEUTROPHIL + PLATELET PATHWAYS
###############################################################################

plate_df <- plate_df %>%
  mutate(
    pathway_wrapped = case_when(
      str_detect(pathway, "RUNX1_REGULATES") ~ str_replace(
        pathway,
        "MEGAKARYOCYTE",
        "\nMEGAKARYOCYTE"  
      ),
      str_detect(pathway, "FACTORS_INVOLVED") ~ str_replace(
        pathway,
        "MEGAKARYOCYTE",
        "\nMEGAKARYOCYTE"
      ),
      TRUE ~ pathway
    )
  )

neuplap <- ggplot(plate_df, aes(x = factor(dataset,levels = tissue_order), y = pathway_wrapped)) +
  geom_point(aes(size = abs(NES), fill = Tissue), 
             shape = 21, stroke = 0.5) +
  theme_bw() + 
  theme(
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = .65),
    axis.text.x = element_text(size = 9, face = "bold",angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9)
    )+
  labs(x="", y = "",title = "Platelets Pathways across SSc datasets", size = "abs(NES)")+
  scale_fill_manual(
    values = tissue_cols,
    name = "Tissue",
    guide = guide_legend(override.aes = list(size = 6))
  ) +
  scale_size_continuous(name = "NES - FDR 10%")

ggsave(plot = neuplap,filename = "results/plots4publication/AllpathsPla_v2.png",
       dpi = 300,width = 8,height = 3.5)



###############################################################################
# DEG ANALYSIS — NEUTROPHIL & PLATELET GENES
###############################################################################

# Gene lists
neutrophil_genes <- c(
  "MPO","ELANE","PADI4","CTSG","CYBB",
  "SELPLG","MMP8","MMP9","HMGB1",
  "S100A8","S100A9","PTX3"
)

platelet_genes <- c(
  "FYN","LYN","SYK","PLCG2",
  "SLP76","FCER1G","ADAM10"
)

# Dataset → tissue mapping
tissue_map <- tibble(
  dataset = c("GSE95065","GSE181549","GSE130955","GSE249550",
              "GSE130953","GSE231691",
              "GSE231693","GSE48149"),
  Tissue = c("Skin","Skin","Skin","Skin",
             "Blood","Blood",
             "Lung","Lung")
)



###############################################################################
# LOAD DEG RESULTS
###############################################################################

dea_res <- list.files("results/dea",
                      pattern = "\\.txt$",
                      full.names = TRUE) %>%
  setNames(gsub("DE_results_|\\.txt$", "",
                basename(.))) %>%
  lapply(read_tsv)


###############################################################################
# FUNCTION — Extract genes from a DEG table
###############################################################################

extract_genes <- function(df, gene_vec, dataset_name){
  df %>%
    filter(SYMBOL %in% gene_vec) %>%
    select(SYMBOL, logFC, adj.P.Val) %>%
    mutate(dataset = dataset_name)
}



###############################################################################
# NEUTROPHIL GENE HEATMAP
###############################################################################

neutro_all <- map2_df(
  dea_res, names(dea_res),
  ~ extract_genes(.x, neutrophil_genes, .y)
)

# Keep max |logFC| if duplicated
neutro_all_clean <- neutro_all %>%
  group_by(SYMBOL, dataset) %>%
  slice_max(abs(logFC), n = 1) %>%
  ungroup() %>%
  left_join(tissue_map, by = "dataset")

# Build logFC matrix
heat_mat <- neutro_all_clean %>%
  select(SYMBOL, dataset, logFC) %>%
  pivot_wider(names_from = dataset, values_from = logFC) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

cols_ordered <- intersect(tissue_order, colnames(heat_mat))
heat_mat <- heat_mat[, cols_ordered, drop = FALSE]

# Build p-value matrix
heat_mat_padj <- neutro_all_clean %>%
  select(SYMBOL, dataset, adj.P.Val) %>%
  pivot_wider(names_from = dataset, values_from = adj.P.Val) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()
heat_mat_padj <- heat_mat_padj[, cols_ordered, drop = FALSE]

# Stars matrix
pheatmap_symbols_neu <- matrix("", nrow = nrow(heat_mat_padj),
                               ncol = ncol(heat_mat_padj))
pheatmap_symbols_neu[heat_mat_padj < 0.05]  <- "*"
pheatmap_symbols_neu[heat_mat_padj < 0.01]  <- "**"
pheatmap_symbols_neu[heat_mat_padj < 0.001] <- "***"

# Tissue annotation
ann_col <- neutro_all_clean %>%
  distinct(dataset, Tissue) %>%
  column_to_rownames("dataset")
ann_col <- ann_col[cols_ordered, , drop = FALSE]

make_italic <- function(x) {
  as.expression(lapply(x, function(y) bquote(italic(.(y)))))
}

lim <- max(abs(heat_mat), na.rm = TRUE)

colors <- colorRampPalette(c("blue", "white", "red"))(100)

breaks <- seq(-lim, lim, length.out = 101)

neuheat <- pheatmap(
  heat_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colors,
  breaks = breaks,          
  annotation_col = ann_col,
  annotation_colors = list(Tissue = tissue_cols),
  labels_row = make_italic(rownames(heat_mat)),
  display_numbers = pheatmap_symbols_neu,
  number_format = "%s",
  fontsize_number = 10,
  main = "Neutrophil-related genes logFC across datasets"
)

ggsave("results/plots4publication/Heatmap_NeuGenes_v2.png",
       plot = neuheat, dpi = 300, width = 7.2, height = 6)



###############################################################################
# PLATELET GENE HEATMAP (Same workflow)
###############################################################################

platelets_all <- map2_df(
  dea_res, names(dea_res),
  ~ extract_genes(.x, platelet_genes, .y)
)

platelets_all_clean <- platelets_all %>%
  group_by(SYMBOL, dataset) %>%
  slice_max(abs(logFC), n = 1) %>%
  ungroup() %>%
  left_join(tissue_map, by = "dataset")

heat_mat <- platelets_all_clean %>%
  select(SYMBOL, dataset, logFC) %>%
  pivot_wider(names_from = dataset, values_from = logFC) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

cols_ordered <- intersect(tissue_order, colnames(heat_mat))
heat_mat <- heat_mat[, cols_ordered, drop = FALSE]

heat_mat_padj <- platelets_all_clean %>%
  select(SYMBOL, dataset, adj.P.Val) %>%
  pivot_wider(names_from = dataset, values_from = adj.P.Val) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()
heat_mat_padj <- heat_mat_padj[, cols_ordered, drop = FALSE]

pheatmap_symbols_neu <- matrix("", nrow = nrow(heat_mat_padj),
                               ncol = ncol(heat_mat_padj))
pheatmap_symbols_neu[heat_mat_padj < 0.05]  <- "*"
pheatmap_symbols_neu[heat_mat_padj < 0.01]  <- "**"
pheatmap_symbols_neu[heat_mat_padj < 0.001] <- "***"

ann_col <- platelets_all_clean %>%
  distinct(dataset, Tissue) %>%
  column_to_rownames("dataset")
ann_col <- ann_col[cols_ordered, , drop = FALSE]

plaheat <- pheatmap(
  heat_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colors,
  breaks = breaks,      
  annotation_col = ann_col,
  annotation_colors = list(Tissue = tissue_cols),
  labels_row = make_italic(rownames(heat_mat)),
  display_numbers = pheatmap_symbols_neu,
  number_format = "%s",
  fontsize_number = 10,
  main = "Platelet-related genes logFC across datasets"
)

ggsave("results/plots4publication/Heatmap_PlateletGenes_v2.png",
       plot = plaheat, dpi = 300, width = 7, height = 6)



###############################################################################
# GSEA ENRICHMENT PANEL (2x4)
###############################################################################

# Build C2 Reactome gene sets
c2_gene_sets <- msigdbr(species = "human", category = "C2")

react <- c2_gene_sets %>%
  filter(gs_subcollection == "CP:REACTOME") %>%
  select(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(all.genes = list(unique(gene_symbol))) %>%
  deframe()


# Create a function for enrichment plot
enrichBy <- function(dataset, geares, dataset_name) {
  
  stats_gsea <- setNames(dataset$logFC, dataset$SYMBOL)
  
  padj_val <- geares %>%
    filter(pathway == "REACTOME_NEUTROPHIL_DEGRANULATION") %>%
    pull(padj) %>%
    signif(3)
  
  plotEnrichment(
    react[["REACTOME_NEUTROPHIL_DEGRANULATION"]],
    stats_gsea
  ) +
    labs(
      title = paste0("Neutrophil Degranulation — ", dataset_name),
      subtitle = paste("Adjusted p-value:", padj_val)
    ) +
    theme_bw()
}

# Generate enrichment plots in desired order
enrich_plots <- imap(dea_res,
                     ~ enrichBy(.x, gsea_res[[.y]], dataset_name = .y))

enrich_plots_ordered <- enrich_plots[tissue_order]

panel_2x4 <- wrap_plots(enrich_plots_ordered, ncol = 2)

ggsave("results/plots4publication/EnrichmentPanel_v2.png",
       plot = panel_2x4, dpi = 300, width = 10, height = 9)
