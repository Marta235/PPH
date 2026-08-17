library(data.table)
library(ggplot2)
library('ramify')
library('psych')
library(pheatmap)
library(dplyr)
library(viridis)
library(ggrastr)

library(wesanderson)

input<-snakemake@input[['data']]
umap<-snakemake@input[['umap']]
plot_out<-snakemake@output[['plot_out']]


cinque<-c('ATOH1','DLL1','GFI1','DEFA5','DEFA6')
dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE)
print(colnames(dato))


cinque<-intersect(cinque, colnames(dato))
print(cinque)
dato<-dato[,cinque]


#set saturation 
sat_max<-9
sat_min<-1
dato[dato>sat_max]<-sat_max
dato[dato<sat_min]<-sat_min


chords<-read.table(file = umap,row.names = 'cell_id',sep=",",header = TRUE)

colori<-rev(rainbow(10))[3:10]
data<-merge(chords,dato,by='row.names')
library(ggplot2)
library(patchwork)
library(ggrastr)
library(scales)

# range globale (senza outlier handling)
lims <- c(sat_min, sat_max)

common_scale <- scale_color_gradientn(
  colours = colori,
  limits  = lims,
  guide   = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barheight = unit(18, "mm"),   # <<< colorbar corta
    barwidth  = unit(3, "mm")     # <<< e stretta
  )
)
xlim_all <- range(data$umap1, na.rm = TRUE)
ylim_all <- range(data$umap2, na.rm = TRUE)

base_theme <- theme_classic(base_size = 7) +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),

    axis.line  = element_line(linewidth = 0.3),

    legend.position = "none",

    aspect.ratio = 1
  )
size_point <- 0.0006

atoh <- ggplot(data, aes(umap1, umap2, color = ATOH1)) +
  rasterize(geom_point(size = size_point), dpi = 300) +
  common_scale +
  base_theme +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = TRUE)

dll<- ggplot(data, aes(umap1, umap2, color = DLL1)) +
  rasterize(geom_point(size = size_point), dpi = 300) +
  common_scale +
  base_theme +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = TRUE)

gfi<-ggplot(data, aes(umap1, umap2, color = GFI1)) +
  rasterize(geom_point(size = size_point), dpi = 300) +
  common_scale +
  base_theme +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = TRUE)

defa5<-ggplot(data, aes(umap1, umap2, color = DEFA5)) +
  rasterize(geom_point(size = size_point), dpi = 300) +
  common_scale +
  base_theme +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = TRUE)

defa6<-ggplot(data, aes(umap1, umap2, color = DEFA6)) +
  rasterize(geom_point(size = size_point), dpi = 300) +
  common_scale +
  base_theme +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = TRUE)


final_plot <- patchwork::wrap_plots(
  atoh, dll, gfi, defa5, defa6,
  ncol = 5
)

ggsave(plot_out, final_plot, width = 180, height = 38, units = "mm", dpi = 300)

