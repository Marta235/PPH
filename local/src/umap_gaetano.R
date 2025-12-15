library(data.table)
library(ggplot2)
library('ramify')
library('psych')
library(pheatmap)
library(dplyr)
library(viridis)

library(wesanderson)

input<-snakemake@input[['data_gae']]
umap<-snakemake@input[['umap']]
plot_out<-snakemake@output[['plot_out']]

dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE)
print(head(dato))
dato$prediction<-as.factor(dato$prediction)

chords<-read.table(file = umap,row.names = 1,sep=",",header = TRUE)
print(head(chords))
colori<-rev(rainbow(10))[3:10]

data<-merge(chords,dato,by='row.names')
print(head(data))
j<-ggplot(data, aes(x=x, y=y,color=prediction)) + 
  geom_point(size=1)+xlab('umap1')+ylab('umap2')+
theme_classic()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())

#fine <- ggarrange( p,k, ncol = 2, common.legend = FALSE)
#ggsave(plot_out, plot=j, width=90, height=90, units="mm")
pdf(plot_out)
print(j)
graphics.off()