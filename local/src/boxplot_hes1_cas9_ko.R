library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)

genes <- c("HES1", "DLL1", "ATOH1")

cas9_path <- snakemake@input$cas9_fpkm
ko_path <- snakemake@input$clones_fpkm
output_plot <- snakemake@output$outplot


# =============================================================================
# CRC0322 Cas9
# =============================================================================

data <- read.table(
  cas9_path,
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1
)

colnames_list <- colnames(data)
split_names <- strsplit(colnames_list, "_")

metadata <- data.frame(
  id = colnames_list,
  sample = sapply(split_names, function(x) x[1]),
  geno = sapply(split_names, function(x) x[2]),
  trattamento = sapply(split_names, function(x) x[3]),
  replica = sapply(split_names, function(x) gsub("\\D", "", x[4])),
  stringsAsFactors = FALSE
)

rownames(metadata) <- metadata$id
metadata$id <- NULL

data <- data[genes, , drop = FALSE]
data <- t(data)

data <- merge(
  data,
  metadata,
  by = "row.names"
)

colnames(data)[1] <- "id"

data$gruppo <- "CRC0322 Cas9"

data$trattamento <- factor(
  data$trattamento,
  levels = c("EGF0.1", "CETUX")
)


# =============================================================================
# CRC0322 ATOH1-KO clones
# =============================================================================

cloni <- read.table(
  ko_path,
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1
)

colnames_list <- colnames(cloni)
split_names <- strsplit(colnames_list, "_")

metacloni <- data.frame(
  id = colnames_list,
  sample = sapply(split_names, function(x) x[1]),
  geno = sapply(split_names, function(x) x[2]),
  trattamento = sapply(split_names, function(x) x[3]),
  replica = sapply(split_names, function(x) gsub("\\D", "", x[4])),
  stringsAsFactors = FALSE
)

rownames(metacloni) <- metacloni$id
metacloni$id <- NULL

cloni <- cloni[genes, , drop = FALSE]
cloni <- t(cloni)

cloni <- merge(
  cloni,
  metacloni,
  by = "row.names"
)

colnames(cloni)[1] <- "id"

cloni$gruppo <- "CRC0322 KO clones"

cloni$trattamento <- factor(
  cloni$trattamento,
  levels = c("EGF0.1", "CETUX")
)


# =============================================================================
# Funzione per definire intervalli regolari dell'asse y
# =============================================================================

nice_breaks_0 <- function(y, k = 5) {
  
  y_max <- max(y, na.rm = TRUE)
  
  if (!is.finite(y_max) || y_max <= 0) {
    return(
      list(
        breaks = c(0, 1),
        top = 1
      )
    )
  }
  
  raw_step <- y_max / k
  mag <- 10^floor(log10(raw_step))
  base <- raw_step / mag
  
  nice_base <- c(1, 2, 2.5, 5, 10)
  step <- nice_base[which(nice_base >= base)[1]] * mag
  
  top <- ceiling(y_max / step) * step
  breaks <- seq(0, top, by = step)
  
  list(
    breaks = breaks,
    top = top
  )
}


# =============================================================================
# Unione dei dati
# =============================================================================

plot_data <- rbind(data, cloni)

plot_data$gruppo <- factor(
  plot_data$gruppo,
  levels = c(
    "CRC0322 Cas9",
    "CRC0322 KO clones"
  )
)

plot_data$trattamento <- factor(
  plot_data$trattamento,
  levels = c(
    "EGF0.1",
    "CETUX"
  )
)

br <- nice_breaks_0(
  plot_data$HES1,
  k = 5
)


# =============================================================================
# Colori associati al trattamento
# =============================================================================

my_colors <- c(
  "EGF0.1" = "black",
  "CETUX" = "red"
)


# =============================================================================
# Plot
# =============================================================================

p_unico <- ggplot(
  plot_data,
  aes(
    x = gruppo,
    y = HES1,
    color = trattamento,
    group = interaction(gruppo, trattamento)
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.55,
    fill = "white",
    position = position_dodge(width = 0.7),
    show.legend = TRUE
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.12,
      jitter.height = 0,
      dodge.width = 0.7
    ),
    size = 0.5,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = my_colors,
    breaks = c("EGF0.1", "CETUX"),
    labels = c(
      "EGF0.1" = "EGF 0.1 ng/ml",
      "CETUX" = "Cetuximab"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "CRC0322 Cas9" = "Cas9",
      "CRC0322 KO clones" = "KO clones"
    ),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_y_continuous(
    limits = c(0, br$top),
    breaks = br$breaks,
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "HES1 (FPKM)",
    color = NULL
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        fill = "white",
        linewidth = 0.6
      )
    )
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 8),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),

    plot.title = element_text(
      hjust = 0.5,
      size = 8,
      margin = margin(b = 0)
    ),

    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box.just = "center",
    legend.key.width = grid::unit(5, "mm"),
    legend.key.height = grid::unit(4, "mm"),
    legend.spacing.x = grid::unit(1, "mm"),
    legend.box.spacing = grid::unit(0, "mm"),
    legend.margin = margin(0, 0, 0, 0),

    axis.text.x = element_text(size = 8),

    plot.margin = margin(
      t = 1,
      r = 2,
      b = 1,
      l = 1,
      unit = "mm"
    )
  )

print(p_unico)

ggsave(
  output_plot,
  p_unico,
  width = 55,
  height = 55,
  units = "mm"
)