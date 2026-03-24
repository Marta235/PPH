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


dato[dato$isPaneth == 'filtered', ]$isPaneth <- 'Others'


dato$prediction <- ifelse(dato$isPaneth == "Paneth", "SPC", "NSPC")
dato$prediction <- factor(dato$prediction, levels = c("SPC", "NSPC"))

chords <- read.table(file = umap, row.names = 'cell_id', sep = ",", header = TRUE)


data <- merge(chords, dato, by = "row.names")

dim_pallini=0.01

j <- ggplot() +
  rasterize(
    geom_point(
      data = data[data$prediction == "NSPC", ],
      aes(x = umap1, y = umap2, color = prediction),
      size = dim_pallini
    ), dpi = 300
  ) +
  rasterize(
    geom_point(
      data = data[data$prediction == "SPC", ],
      aes(x = umap1, y = umap2, color = prediction),
      size = dim_pallini
    ), dpi = 300
  ) +
  scale_color_manual(
    name = '',
    values = c(
      "SPC" = "#099963",
      "NSPC"   = "#76069A"
    )
  ) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 1)
  ))

  

ggsave(plot_out, plot = j, width = 75, height = 58, units = "mm")
