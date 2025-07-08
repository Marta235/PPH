library(data.table)

input<-snakemake@input[['data']]

out<-snakemake@output[['out']]
#leggo dato
dato<- read.table(file = input,row.names = 1,sep=",",header = TRUE)
dato_t<- transpose(dato)
rownames(dato_t) <- colnames(dato)
colnames(dato_t)<-rownames(dato)

stringa_split<-function(stringa){
  res<-strsplit(stringa , split = ":")[[1]][2]
  return(res)
}
colnames(dato_t) <- sapply(colnames(dato_t),FUN=stringa_split)
#cinque<-c()
#cinque<-c('TSPAN6','TNMD','GCLC','ENPP4','SEMA3F')
cinque<-c("ATOH1","LEF1","GFI1","DLL1","DEFA5","DEFA6","OLFM4","NOTCH1","NOTCH2","HES1","DLL4","APCDD1","WNT6","LGR5","BTC","AREG","EGF","ERBB2","ERBB3","EREG","HBEGF","TGFA","EGFR","ERBB2","SPDEF","CREB3L4")
#cinque<-read.table(file = gene_list,sep=",",header = FALSE)
#cinque<-cinque$V1
#cinque<-intersect(cinque, colnames(dato_t))

cinque_df<-dato_t[,cinque]

write.table(cinque_df,file=out,sep=',',quote=FALSE)
