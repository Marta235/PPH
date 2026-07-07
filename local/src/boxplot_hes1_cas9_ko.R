library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
genes <- c('HES1','DLL1','ATOH1')

cas9_path<-snakemake@input$cas9_fpkm
ko_path<-snakemake@input$clones_fpkm
output_plot<-snakemake@output$outplot

data <- read.table(cas9_path, stringsAsFactors = FALSE, header = TRUE, row.names = 1)

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
data <- merge(data, metadata, by = "row.names")
colnames(data)[1] <- "id"

data$gruppo <- "CRC0322 Cas9"
data$trattamento <- factor(data$trattamento, levels = c("EGF0.1", "CETUX"))

cloni <- read.table(ko_path, stringsAsFactors = FALSE, header = TRUE, row.names = 1)

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
cloni <- merge(cloni, metacloni, by = "row.names")
colnames(cloni)[1] <- "id"

cloni$gruppo <- "CRC0322 KO clones"
cloni$trattamento <- factor(cloni$trattamento, levels = c("EGF0.1", "CETUX"))


nice_breaks_0 <- function(y, k = 5) {
  y_max <- max(y, na.rm = TRUE)
  
  if (!is.finite(y_max) || y_max <= 0) {
    return(list(breaks = c(0, 1), top = 1))
  }
  
  raw_step <- y_max / k
  mag <- 10^floor(log10(raw_step))
  base <- raw_step / mag
  
  nice_base <- c(1, 2, 2.5, 5, 10)
  step <- nice_base[which(nice_base >= base)[1]] * mag
  
  top <- ceiling(y_max / step) * step
  breaks <- seq(0, top, by = step)
  
  list(breaks = breaks, top = top)
}

plot_data <- rbind(data, cloni)
plot_data$gruppo <- factor(plot_data$gruppo, levels = c("CRC0322 Cas9", "CRC0322 KO clones"))

br <- nice_breaks_0(plot_data$HES1, k = 5)
my_colors <- c(
  "CRC0322 Cas9"      = "black",
  "CRC0322 KO clones" = "#3F6FA8"
)
p_unico <- ggplot(plot_data, aes(x = trattamento, y = HES1, color = gruppo)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7,
    position = position_dodge(width = 0.8),
    alpha = 0.35
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.18,
      jitter.height = 0,
      dodge.width = 0.8
    ),
    size = 0.3) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors)+
  theme_bw() +
  theme(
    text = element_text(size = 8),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.title = element_blank()
  ) +
  scale_x_discrete(labels = c("EGF0.1" = "EGF 0.1 ng/ml", "CETUX" = "Cetuximab")) +
  scale_y_continuous(
    limits = c(0, br$top),
    breaks = br$breaks,
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "HES1 (FPKM)")

print(p_unico)
ggsave(output_plot, p_unico, width = 90, height = 55, units = "mm")

