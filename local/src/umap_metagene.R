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


cinque<-intersect(cinque, colnames(dato))
print(cinque)
dato<-dato[,cinque]


metagene_cinque<-apply(dato,1,mean)
metagene_cinque<-as.data.frame(metagene_cinque)
rownames(metagene_cinque)<-rownames(dato)
colnames(metagene_cinque)<-'metagene'
#set saturation to 10
sat_max<-as.numeric(snakemake@wildcards[['sat_max']])
sat_min<-as.numeric(snakemake@wildcards[['sat_min']])
metagene_cinque[metagene_cinque>sat_max]<-sat_max
metagene_cinque[metagene_cinque<sat_min]<-sat_min
#write.csv(metagene_cinque, meta_ou, row.names=TRUE)

chords<-read.table(file = umap,row.names = 'cell_id',sep=",",header = TRUE)
# print(head(metagene_cinque))
# print(rownames(chords))
colori<-rev(rainbow(10))[3:10]
data<-merge(chords,metagene_cinque,by='row.names')

xr <- range(data$umap1, na.rm = TRUE)
yr <- range(data$umap2, na.rm = TRUE)
r <- max(diff(xr), diff(yr)) / 2
cx <- mean(xr)
cy <- mean(yr)

j <- ggplot(data, aes(x = umap1, y = umap2, color = metagene)) +
  rasterize(geom_point(size = 0.01), dpi = 300) +
  scale_color_gradientn(colours = colori) +
  labs(color = "SPC metagene log2cpm") +
  xlab("UMAP1") + ylab("UMAP2") +
  coord_fixed(
    xlim = c(cx - r, cx + r),
    ylim = c(cy - r, cy + r)
  ) +
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
  guides(color = guide_colorbar(
    barheight = unit(11, "mm"),
    barwidth  = unit(2, "mm")
  ))
#fine <- ggarrange( p,k, ncol = 2, common.legend = FALSE)
ggsave(plot_out, plot = j, width = 70, height = 58, units = "mm")
