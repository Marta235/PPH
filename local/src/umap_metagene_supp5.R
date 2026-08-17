library(data.table)
library(ggplot2)
library('ramify')
library('psych')
library(pheatmap)
library(dplyr)
library(viridis)
library(ggrastr)
library(patchwork)
library(scales)

library(wesanderson)

input_nt<-snakemake@input[['data_nt']]
umap_nt<-snakemake@input[['umap_nt']]
plot_out<-snakemake@output[['plot_out']]

input_cet<-snakemake@input[['data_cet']]
umap_cet<-snakemake@input[['umap_cet']]

cinque<-c('SMAD1','FGFR3')
dato_nt<- read.table(file = input_nt,row.names = 1,sep=",",header = TRUE)
cinque<-intersect(cinque, colnames(dato_nt))
dato_nt<-dato_nt[,cinque]

dato_cet<- read.table(file = input_cet,row.names = 1,sep=",",header = TRUE)
cinque<-intersect(cinque, colnames(dato_cet))
dato_cet<-dato_cet[,cinque]
  




#set saturation to 10
sat_max<-as.numeric(snakemake@wildcards[['sat_max']])
sat_min<-as.numeric(snakemake@wildcards[['sat_min']])
dato_nt[dato_nt>sat_max]<-sat_max
dato_nt[dato_nt<sat_min]<-sat_min
dato_cet[dato_cet>sat_max]<-sat_max
dato_cet[dato_cet<sat_min]<-sat_min

#write.csv(metagene_cinque, meta_ou, row.names=TRUE)

chords_nt<-read.table(file = umap_nt,row.names = 'cell_id',sep=",",header = TRUE)
chords_cet<-read.table(file = umap_cet,row.names = 'cell_id',sep=",",header = TRUE)

colori<-rev(rainbow(10))[3:10]
data_nt<-merge(chords_nt,dato_nt,by='row.names')
data_cet<-merge(chords_cet,dato_cet,by='row.names')

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
xlim_all <- range(c(data_nt$umap1, data_cet$umap1), na.rm = TRUE)
ylim_all <- range(c(data_nt$umap2, data_cet$umap2), na.rm = TRUE)

base_theme <- theme_classic(base_size = 7) +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 7),
    axis.line  = element_line(linewidth = 0.3),
    plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"),
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 6),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    aspect.ratio = 1
  )

smad_nt <- ggplot(data_nt, aes(umap1, umap2, color = SMAD1)) +
  rasterize(geom_point(size = 0.08), dpi = 300) +
  common_scale +
  labs( x = "umap1", y = "umap2", color = "Expression") +#title = "SMAD1",
  base_theme +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = TRUE)

smad_cet <- ggplot(data_cet, aes(umap1, umap2, color = SMAD1)) +
  rasterize(geom_point(size = 0.08), dpi = 300) +
  common_scale +
  labs( x = "umap1", y = "umap2", color = "Expression") + #title = "SMAD1",
  base_theme +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = TRUE)

fgfr3_nt <- ggplot(data_nt, aes(umap1, umap2, color = FGFR3)) +
  rasterize(geom_point(size = 0.08), dpi = 300) +
  common_scale +
  labs( x = "umap1", y = "umap2", color = "Expression") +#title = "FGFR3",
  base_theme +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = TRUE)

fgfr3_cet <- ggplot(data_cet, aes(umap1, umap2, color = FGFR3)) +
  rasterize(geom_point(size = 0.08), dpi = 300) +
  common_scale +
  labs( x = "umap1", y = "umap2", color = "Expression") + #title = "FGFR3",
  base_theme +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = TRUE)


final_plot <- patchwork::wrap_plots(
  smad_nt, smad_cet, fgfr3_nt, fgfr3_cet,
  ncol = 4,
  guides = "collect"
)

ggsave(plot_out, final_plot, width = 170, height = 42, units = "mm", dpi = 300)

