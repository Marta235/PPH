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
#dato <- dato[dato$isPaneth != 'filtered', ]


dato$prediction<-dato$isPaneth
dato$prediction <- factor(dato$prediction, levels = c("Paneth", "Others","filtered"))

chords <- read.table(file = umap, row.names = 'cell_id', sep = ",", header = TRUE)


data <- merge(chords, dato, by = "row.names")

dim_pallini=0.1
# 4) Plot
j <- ggplot() +
  rasterize(
    geom_point(
      data = data[data$prediction == "Others", ],
      aes(x = umap1, y = umap2, color = prediction),
      size = dim_pallini
    ), dpi = 300
  ) +
  rasterize(
    geom_point(
      data = data[data$prediction == "Paneth", ],
      aes(x = umap1, y = umap2, color = prediction),
      size = dim_pallini
    ), dpi = 300
  ) +
  rasterize(
    geom_point(
      data = data[data$prediction == "filtered", ],
      aes(x = umap1, y = umap2, color = prediction),
      size = dim_pallini
    ), dpi = 300
  ) +
  scale_color_manual(
    name = '',
    values = c(
      "Paneth" = "#099963",
      "Others"   = "#76069A",
      "filtered" = "goldenrod"
    )
  ) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_blank(),
    #legend.position = "none"
  )

  

ggsave(plot_out, plot = j, width = 120, height = 100, units = "mm")
