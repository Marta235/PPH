
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
kmeans_out<-snakemake@output[['out']]
log_f <- snakemake@log[['log']]
sink(log_f)

dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE)
dato_t<- transpose(dato)
rownames(dato_t) <- colnames(dato)
colnames(dato_t)<-rownames(dato)

stringa_split<-function(stringa){
  res<-strsplit(stringa , split = ":")[[1]][2]
  return(res)
}

colnames(dato_t) <- sapply(colnames(dato_t),FUN=stringa_split)
# Select Paneth cell markers
cinque<-c("ATOH1","GFI1","DLL1","DEFA5","DEFA6")
cinque_df<-dato_t[,cinque]
# compute cluster with kmeans
cl <- kmeans(cinque_df, 2)

meta_mu<-apply(cl$centers,1,geometric.mean)
ordine<-order(unlist(meta_mu))
ordinate<-c('Others','Paneth')

final_ordinate<-c()
for (i in seq(1,2)){
  final_ordinate[ordine[i]] <-ordinate[i] 
}
cluster_id<-cl$cluster

isPaneth<-c()
for (el in cluster_id){
  isPaneth<-append(isPaneth,final_ordinate[el])
}
posteriors<-cbind(cluster_id,isPaneth)

posteriors<-data.frame(posteriors)
posteriors$cluster_id <- as.numeric(posteriors$cluster_id) 

conteggio <- table(posteriors$isPaneth)
print(conteggio)

sink()
write.table(posteriors,file=kmeans_out,sep=',',quote=FALSE)



