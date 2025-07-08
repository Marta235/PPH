library(data.table)
library(mclust)
library(cluster)
library(reshape2)
library(ggplot2)
library('ramify')
library('psych')
library(pheatmap)
library(resample)
library(dplyr)
library(viridis)
library(igraph)
library(patchwork)

input<-snakemake@input[['data']]
silu_out<-snakemake@output[['out']]
log_f <- snakemake@log[['log']]
clust_path<-snakemake@input[['clusters']]
sink(log_f)
posteriors <- read.table(clust_path,sep=',',header=TRUE,row.names=1)
dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE)
dato_t<- transpose(dato)
rownames(dato_t) <- colnames(dato)
colnames(dato_t)<-rownames(dato)
stringa_split<-function(stringa){
  res<-strsplit(stringa , split = ":")[[1]][2]
  return(res)
}
colnames(dato_t) <- sapply(colnames(dato_t),FUN=stringa_split)
library(vegan)
posteriors<-posteriors[rownames(dato_t),]
sil<-silhouette(posteriors$cluster_id,dist(dato_t))

silu<-as.data.frame.matrix(sil)
rownames(silu)<-rownames(posteriors)

silu$preSilh<-posteriors$isPaneth

silu$postSilh<-posteriors$isPaneth
silu$preSilh<-as.character(silu$preSilh)
silu$postSilh<-as.character(silu$postSilh)
silu$postSilh[silu$sil_width<0 & silu$preSilh=='Paneth']<-'filtered'

silu$isPaneth<-silu$postSilh
silu$postSilh[silu$postSilh=='filtered']<-0
silu$postSilh[silu$postSilh=='Others']<-0
silu$postSilh[silu$postSilh=='Paneth']<-1

silu$isPaneth<-as.factor(silu$isPaneth)

conteggio <- table(silu$isPaneth)
print(conteggio)

sink()
write.table(silu,file=silu_out,sep=',',quote=FALSE)

