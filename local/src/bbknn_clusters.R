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
cluster<-snakemake@input[['cluster']]
plot_out<-snakemake@output[['plot_out']]

dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
clusters<- read.table(file = cluster,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
cluster<-cluster[cluster$sample==snakemake@wildcards[['sample']],]
umap_data<- read.table(file = umap,row.names = 'uid',sep=",",header = TRUE,stringsAsFactors = FALSE)

to_plot=c(colnames(dato),'leiden_r0.3','phase')
#set saturation to 10
#sat_max<-as.numeric(snakemake@wildcards[['sat_max']])
#sat_min<-as.numeric(snakemake@wildcards[['sat_min']])
#metagene_cinque[metagene_cinque>sat_max]<-sat_max
#metagene_cinque[metagene_cinque<sat_min]<-sat_min
#write.csv(metagene_cinque, meta_ou, row.names=TRUE)


colori<-rev(rainbow(10))[3:10]
data<-merge(umap_data,dato,by='row.names')
print(head(data))
data<-merge(data,clusters,by.x='Row.names',by.y='cell_id')
pdf(plot_out)
for (col in to_plot)){
j<-ggplot(data, aes(x=umap1, y=umap2,color=data[[col]])) + 
rasterize(geom_point(size=1),dpi=300)+scale_color_gradientn(colours = colori)+
labs(title = gene, color=gene)+xlab('umap1')+ylab('umap2')+
theme_classic()+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())
print(j)
}
graphics.off()