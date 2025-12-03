library(data.table)
library(ggplot2)
library(ramify)
library(psych)
library(pheatmap)
library(dplyr)
library(viridis)
library(ggrastr)
library(wesanderson)

input    <- snakemake@input[['data_kmeans']]
umap     <- snakemake@input[['umap']]
plot_out <- snakemake@output[['plot_out']]

dato <- read.table(file = input, row.names = 1, sep = ",", header = TRUE)
print(head(dato))

# 1) Togliamo SUBITO le celle 'filtered'
dato <- dato[dato$isPaneth != 'filtered', ]

# 2) Creiamo direttamente una colonna prediction pulita
#    Paneth -> PCL cell, tutto il resto -> Others
dato$prediction <- ifelse(dato$isPaneth == "Paneth", "PCL cell", "Others")
dato$prediction <- factor(dato$prediction, levels = c("PCL cell", "Others"))

chords <- read.table(file = umap, row.names = 'cell_id', sep = ",", header = TRUE)

# 3) Un solo merge è sufficiente
data <- merge(chords, dato, by = "row.names")

# Giusto per controllo:
print(unique(data$prediction))
print(levels(data$prediction))  # qui NON ci deve essere 'filtered'

# 4) Plot
j <- ggplot(data, aes(x = umap1, y = umap2, color = prediction)) + 
  rasterize(geom_point(size = 0.1), dpi = 300) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  scale_color_manual(
    name = '',
    values = c(
      "PCL cell" = "#1D7937",
      "Others"   = "#772B84"
    )
  ) +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_blank(),legend.position = "none")
  

ggsave(plot_out, plot = j, width = 100, height = 100, units = "mm")
