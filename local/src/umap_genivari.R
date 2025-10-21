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

dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)

print(colnames(dato))

#set saturation to 10
#sat_max<-as.numeric(snakemake@wildcards[['sat_max']])
#sat_min<-as.numeric(snakemake@wildcards[['sat_min']])
#metagene_cinque[metagene_cinque>sat_max]<-sat_max
#metagene_cinque[metagene_cinque<sat_min]<-sat_min
#write.csv(metagene_cinque, meta_ou, row.names=TRUE)

chords<-read.table(file = umap,row.names = 1,sep=",",header = TRUE)
colori<-rev(rainbow(10))[3:10]
data<-merge(chords,dato,by='row.names')

pdf(plot_out)
for (gene in colnames(dato)){
j<-ggplot(data, aes(x=x, y=y,color=data[[gene]])) + 
geom_point(size=1)+scale_color_gradientn(colours = colori)+
labs(title = gene, color=gene)+xlab('umap1')+ylab('umap2')+
theme_classic()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())
print(j)
}
graphics.off()