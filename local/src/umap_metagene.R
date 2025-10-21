library(data.table)
library(ggplot2)
library('ramify')
library('psych')
library(pheatmap)
library(dplyr)
library(viridis)

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

chords<-read.table(file = umap,row.names = 1,sep=",",header = TRUE)
colori<-rev(rainbow(10))[3:10]
data<-merge(chords,metagene_cinque,by='row.names')
j<-ggplot(data, aes(x=x, y=y,color=metagene)) + 
  geom_point(size=1)+scale_color_gradientn(colours = colori,limits=c(sat_min,sat_max))+#+scale_color_viridis(limits=c(0,10),direction = -1)
labs(title = "metagene", color="metagene")+xlab('umap1')+ylab('umap2')+
theme_classic()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())

#fine <- ggarrange( p,k, ncol = 2, common.legend = FALSE)
ggsave(plot_out, plot=j, width=90, height=90, units="mm")
pdf(plot_out)
print(j)
graphics.off()