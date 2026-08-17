library(ComplexHeatmap)
library(circlize)
library(grid)

gsea_path <- snakemake@input[["gsea_cas9"]]
gsea_path_cloni <- snakemake@input[["gsea_cloni"]]
plot_out <- snakemake@output[["plot_out"]]


########################################
# GSEA ANALYSIS NOTCH
########################################

data_gsea <- read.table(
  gsea_path,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

data_gsea_cloni <- read.table(
  gsea_path_cloni,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)


# NOTCH pathways enriched positively in CAS9
notch_gsea <- data_gsea[
  grep("NOTCH", data_gsea$ID),
]

notch_gsea <- notch_gsea[
  notch_gsea$enrichmentScore > 0,
]


########################################
# Build table matching pathways by ID
########################################

cas9 <- notch_gsea[, c(
  "ID",
  "NES",
  "pvalue",
  "p.adjust"
)]

colnames(cas9) <- c(
  "ID",
  "NES_cas9",
  "pval_cas9",
  "p.adjust_cas9"
)


cloni <- data_gsea_cloni[, c(
  "ID",
  "NES",
  "pvalue",
  "p.adjust"
)]

colnames(cloni) <- c(
  "ID",
  "NES_cloni",
  "pval_cloni",
  "p.adjust_cloni"
)


# Manteniamo solo le pathways NOTCH selezionate in CAS9
per_heatmap <- merge(
  cas9,
  cloni,
  by = "ID",
  all.x = TRUE,
  sort = FALSE
)

rownames(per_heatmap) <- per_heatmap$ID
per_heatmap$ID <- NULL


########################################
# NES heatmap
########################################

nes_cols <- c(
  "NES_cas9",
  "NES_cloni"
)

mat_nes <- as.matrix(
  per_heatmap[, nes_cols, drop = FALSE]
)


# order based on CAS9 NES
ord <- order(
  mat_nes[, "NES_cas9"],
  decreasing = TRUE
)


col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#4575b4", "white", "#d73027")
)


ht_nes <- Heatmap(
  mat_nes,
  name = "NES",
  col = col_fun,

  row_order = ord,

  cluster_columns = FALSE,
  cluster_rows = FALSE,

  show_row_dend = FALSE,
  show_row_names = TRUE,

  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  row_names_max_width = unit(12, "cm"),

  heatmap_legend_param = list(
    at = c(-2, 0, 2),
    labels = c("-2", "0", "+2")
  )
)


########################################
# adjusted p-value heatmap
########################################

pval_cols <- c(
  "p.adjust_cas9",
  "p.adjust_cloni"
)

mat_pval <- as.matrix(
  per_heatmap[, pval_cols, drop = FALSE]
)


pval_th <- 0.1

mat_pval_cat <- ifelse(
  mat_pval < pval_th,
  "significativo",
  "non significativo"
)

mat_pval_cat[is.na(mat_pval)] <- "NA"

mat_pval_cat <- as.matrix(mat_pval_cat)


colori_pval <- c(
  "significativo" = "goldenrod",
  "non significativo" = "lightgrey",
  "NA" = "white"
)


ht_pval <- Heatmap(
  mat_pval_cat,

  name = "adjusted p-value",

  row_order = ord,

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  show_row_names = FALSE,

  col = colori_pval,

  rect_gp = gpar(
    col = "white",
    lwd = 1
  ),

  heatmap_legend_param = list(
    at = c(
      "significativo",
      "non significativo",
      "NA"
    ),
    labels = c(
      paste0("< ", pval_th),
      paste0(">= ", pval_th),
      "NA"
    )
  )
)


########################################
# SAVE
########################################

pdf(
  plot_out,
  width = 10,
  height = 6
)

draw(
  ht_nes + ht_pval,
  heatmap_legend_side = "right"
)

dev.off()