library(ggplot2)

input <- snakemake@input[["deseq"]]
plot_out <- snakemake@output[["plot"]]

########################################
# READ DATA
########################################

cas9 <- read.table(
  input,
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1
)

cas9$gene <- rownames(cas9)


########################################
# DEFINE SIGNIFICANCE
########################################

lfc_th <- 0.58
padj_th <- 0.05

cas9$significant <- "Downregulated or NS"

cas9$significant[
  cas9$padj < padj_th &
  cas9$log2FoldChange > lfc_th
] <- "Upregulated"

cas9$significant[
  cas9$padj < padj_th &
  cas9$log2FoldChange < -lfc_th
] <- "Downregulated or NS"


########################################
# HANDLE PADJ = 0
########################################

cas9$padj_plot <- cas9$padj

cas9$padj_plot[
  cas9$padj_plot == 0
] <- 1e-235


########################################
# VOLCANO
########################################

volcano <- ggplot(
  cas9,
  aes(
    x = log2FoldChange,
    y = -log10(padj_plot),
    color = significant
  )
) +
  geom_point(size = 1) +

  scale_color_manual(
    values = c(
      "Downregulated or NS" = "#d3d3d3",
      "Upregulated" = "#CF4030"
    )
  ) +

  geom_vline(
    xintercept = c(-lfc_th, lfc_th),
    linetype = "dashed",
    color = "black"
  ) +

  geom_hline(
    yintercept = -log10(padj_th),
    linetype = "dashed",
    color = "black"
  ) +

  scale_x_continuous(
    limits = c(-8, 8),
    breaks = c(-8, -4, 0, 4, 8)
  ) +

  scale_y_continuous(
    limits = c(0, 240),
    breaks = c(0, 80, 180, 240)
  ) +

  theme_minimal() +

  labs(
    x = "log2FC",
    y = "-log10 adjusted p-value",
    color = NULL
  ) +

  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top")
  )


########################################
# SAVE
########################################

ggsave(
  filename = plot_out,
  plot = volcano,
  width = 12,
  height = 10,
  units = "cm",
  dpi = 300
)