library(data.table)
library(ggplot2)
library(uwot)

input<-snakemake@input[['data']]
data_out<-snakemake@output[['out']]
plot_out<-snakemake@output[['umap_out']]

min_dist<-0.5
n_neighbors<-50

dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE)
dato_t<- transpose(dato)
rownames(dato_t) <- colnames(dato)
colnames(dato_t)<-rownames(dato)

umap_result <- umap(dato_t,min_dist = min_dist,n_neighbors = n_neighbors)
umap_result<-as.data.frame(umap_result)
colnames(umap_result)<-c('x','y')
write.csv(umap_result, data_out, row.names = TRUE)
