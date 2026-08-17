library(ggplot2)

cas9_path  <- snakemake@input[["cas9"]]
combo_path <- snakemake@input[["combo"]]
ko_path    <- snakemake@input[["ko"]]

plot_out   <- snakemake@output[["plot"]]
log_out    <- snakemake@log[["fisher"]]


########################################
# READ DATA
########################################

cas9 <- read.table(
  cas9_path,
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1
)

creni <- read.table(
  combo_path,
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1
)

ko <- read.table(
  ko_path,
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1
)


########################################
# THRESHOLDS
########################################

th_padj <- 0.05
th_log2fc <- 0.58


########################################
# UP / DOWN SETS
########################################

cas9_up <- rownames(
  cas9[
    !is.na(cas9$padj) &
    cas9$padj < th_padj &
    cas9$log2FoldChange > th_log2fc,
  ]
)

cas9_down <- rownames(
  cas9[
    !is.na(cas9$padj) &
    cas9$padj < th_padj &
    cas9$log2FoldChange < -th_log2fc,
  ]
)


combo_up <- rownames(
  creni[
    !is.na(creni$padj) &
    creni$padj < th_padj &
    creni$log2FoldChange > th_log2fc,
  ]
)

combo_down <- rownames(
  creni[
    !is.na(creni$padj) &
    creni$padj < th_padj &
    creni$log2FoldChange < -th_log2fc,
  ]
)


ko_up <- rownames(
  ko[
    !is.na(ko$padj) &
    ko$padj < th_padj &
    ko$log2FoldChange > th_log2fc,
  ]
)

ko_down <- rownames(
  ko[
    !is.na(ko$padj) &
    ko$padj < th_padj &
    ko$log2FoldChange < -th_log2fc,
  ]
)


########################################
# UP CONTINGENCY TABLE
########################################

# CAS9-up genes for which data exist in all comparisons
universe_up <- intersect(
  cas9_up,
  intersect(rownames(ko), rownames(creni))
)

# CAS9 genes that are NOT significantly up in combo / KO
cas9_combo_up <- setdiff(universe_up, combo_up)
cas9_ko_up <- setdiff(universe_up, ko_up)


# not in combo AND not in KO
a_up <- intersect(
  cas9_combo_up,
  cas9_ko_up
)

# not in combo, but in KO
b_up <- setdiff(
  cas9_combo_up,
  a_up
)

# in combo, but not in KO
c_up <- setdiff(
  cas9_ko_up,
  a_up
)

# in both combo and KO
d_up <- setdiff(
  universe_up,
  union(
    union(a_up, b_up),
    c_up
  )
)


contingency_table_up <- matrix(
  c(
    length(a_up),
    length(b_up),
    length(c_up),
    length(d_up)
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    "cas9-combo" = c("not_in", "in"),
    "cas9-ko" = c("not_in", "in")
  )
)

fisher_up <- fisher.test(
  contingency_table_up,
  alternative = "greater"
)


########################################
# DOWN CONTINGENCY TABLE
########################################

universe_down <- intersect(
  cas9_down,
  intersect(rownames(ko), rownames(creni))
)

cas9_combo_down <- setdiff(
  universe_down,
  combo_down
)

cas9_ko_down <- setdiff(
  universe_down,
  ko_down
)


a_down <- intersect(
  cas9_combo_down,
  cas9_ko_down
)

b_down <- setdiff(
  cas9_combo_down,
  a_down
)

c_down <- setdiff(
  cas9_ko_down,
  a_down
)

d_down <- setdiff(
  universe_down,
  union(
    union(a_down, b_down),
    c_down
  )
)


contingency_table_down <- matrix(
  c(
    length(a_down),
    length(b_down),
    length(c_down),
    length(d_down)
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    "cas9-combo" = c("not_in", "in"),
    "cas9-ko" = c("not_in", "in")
  )
)

fisher_down <- fisher.test(
  contingency_table_down,
  alternative = "greater"
)


########################################
# TOTAL FISHER
########################################

contingency_table <- matrix(
  c(
    length(a_up) + length(a_down),
    length(b_up) + length(b_down),
    length(c_up) + length(c_down),
    length(d_up) + length(d_down)
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    "cas9-combo" = c("not_in", "in"),
    "cas9-ko" = c("not_in", "in")
  )
)

fisher_total <- fisher.test(
  contingency_table,
  alternative = "greater"
)


########################################
# WRITE FISHER RESULTS TO LOG
########################################

sink(log_out)

cat("Thresholds\n")
cat("padj <", th_padj, "\n")
cat("|log2FC| >", th_log2fc, "\n\n")


cat("====================================\n")
cat("UPREGULATED GENES\n")
cat("====================================\n\n")

print(contingency_table_up)

cat("\nFisher test:\n")
print(fisher_up)

cat("\n\n")


cat("====================================\n")
cat("DOWNREGULATED GENES\n")
cat("====================================\n\n")

print(contingency_table_down)

cat("\nFisher test:\n")
print(fisher_down)

cat("\n\n")


cat("====================================\n")
cat("UP + DOWN\n")
cat("====================================\n\n")

print(contingency_table)

cat("\nFisher test:\n")
print(fisher_total)

sink()


########################################
# SCATTER PLOT
########################################

creni_plot <- creni[, "log2FoldChange", drop = FALSE]
ko_plot <- ko[, "log2FoldChange", drop = FALSE]

colnames(creni_plot) <- "log2FC_ComboVsEGF"
colnames(ko_plot) <- "log2FC_CetuxvsEGF_cloni"


merged <- merge(
  creni_plot,
  ko_plot,
  by = "row.names"
)

rownames(merged) <- merged$Row.names
merged$Row.names <- NULL


# Come nel notebook: il plot considera i geni UP in CAS9
merged <- merged[universe_up, , drop = FALSE]


merged$tipo <- ifelse(
  rownames(merged) %in% a_up,
  "not_in_both_up",
  ifelse(
    rownames(merged) %in% b_up,
    "not_in_combo",
    ifelse(
      rownames(merged) %in% c_up,
      "not_in_KO",
      "as_cas9"
    )
  )
)


merged$tipo <- factor(
  merged$tipo,
  levels = c(
    "as_cas9",
    "not_in_combo",
    "not_in_KO",
    "not_in_both_up"
  )
)

merged <- merged[
  order(as.numeric(merged$tipo)),
  ,
  drop = FALSE
]


scatter_plot <- ggplot(
  merged,
  aes(
    x = log2FC_CetuxvsEGF_cloni,
    y = log2FC_ComboVsEGF,
    color = tipo,
    alpha = tipo == "as_cas9"
  )
) +

  geom_point(
    size = 1,
    stroke = 0
  ) +

  labs(
    x = "log2FC Cetux vs EGF KO clones",
    y = "log2FC Combo vs EGF"
  ) +

  theme_minimal() +

  scale_color_manual(
    values = c(
      "not_in_both_up" = "#9CAFB7",
      "not_in_combo"   = "#254E70",
      "not_in_KO"      = "#5D2A42",
      "as_cas9"        = "#CF4030"
    )
  ) +

  scale_alpha_manual(
    values = c(
      "TRUE" = 1,
      "FALSE" = 1
    ),
    guide = "none"
  ) +

  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.5
  ) +

  geom_vline(
    xintercept = 0,
    color = "black",
    linewidth = 0.5
  ) +

  geom_abline(
    intercept = 0,
    slope = 1,
    color = "black",
    linetype = "dashed",
    linewidth = 0.5
  ) +

  geom_abline(
    intercept = 0,
    slope = -1,
    color = "black",
    linetype = "dashed",
    linewidth = 0.5
  ) +

  geom_hline(
    yintercept = c(-th_log2fc, th_log2fc),
    color = "black",
    linetype = "dashed",
    linewidth = 0.5
  ) +

  geom_vline(
    xintercept = c(-th_log2fc, th_log2fc),
    color = "black",
    linetype = "dashed",
    linewidth = 0.5
  ) +

  scale_x_continuous(
    limits = c(-2.5, 5),
    breaks = c(-2.5, 0, 2.5, 5)
  ) +

  scale_y_continuous(
    limits = c(-2.5, 5),
    breaks = c(-2.5, 0, 2.5, 5)
  ) +

  theme(
    panel.grid = element_blank(),

    axis.line = element_line(
      color = "black"
    ),

    axis.ticks = element_line(
      color = "black"
    ),

    legend.position = c(0.98, 0.02),

    legend.justification = c(
      "right",
      "bottom"
    ),

    legend.background = element_rect(
      fill = "white",
      color = NA
    ),

    legend.key = element_blank()
  )


########################################
# SAVE
########################################

ggsave(
  filename = plot_out,
  plot = scatter_plot,
  width = 17,
  height = 16,
  units = "cm",
  dpi = 300
)