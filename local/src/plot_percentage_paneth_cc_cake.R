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
library(wesanderson)

library(reshape2) 
library(reshape) 
library(tidyr)


clust<-snakemake@input[['clust']]
cc<-snakemake@input[['cc']]
plot_out<-snakemake@output[['plot_out']]
log_out<-snakemake@log[['log']]


dato<- read.table(file = clust,row.names = 1,sep=",",header = TRUE)


dato$isPaneth<-as.character(dato$isPaneth)
dato$isPaneth[dato$isPaneth=='filtered']<-'Others'
dato$isPaneth[dato$isPaneth=='nPaneth']<-'Others'

dato$isPaneth<-as.factor(dato$isPaneth)



dato_cc<- read.table(file = cc,row.names = 1,sep=",",header = FALSE)
dato<-merge(dato,dato_cc,by='row.names')

dato.summary <- dato %>% group_by(isPaneth) %>% 
  summarise(total_count=n(),.groups = 'drop') %>% 

  mutate(percent =total_count/sum(total_count))


dato.summary_paneth <-dato %>% group_by(isPaneth,V2) %>% 
  summarise(total_count=n(),.groups = 'drop') %>% 
  group_by(isPaneth) %>% 
  mutate(percent =total_count/sum(total_count),mucca=1-( cumsum(percent) - 0.5*percent))
colori<-c("#76069A","#099963")
paneth <- ggplot(dato.summary, aes(x = isPaneth, y = percent, fill = isPaneth)) +
  geom_bar(stat = "identity", width = .7, fill = colori, lwd = 0.1) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(y = "", x = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),   # rimuove griglia maggiore
    panel.grid.minor = element_blank(),   # rimuove griglia minore
    panel.border     = element_blank(),   # rimuove bordo del pannello
    axis.line        = element_line(color = "black", size = 0.5), # mostra assi
    legend.position  = "none"
  )

#cet
cont<- dato.summary_paneth[,c('isPaneth','V2','total_count')]
spread_df <- spread(cont, key = V2, value = total_count, fill = 0)
spread_df<-as.data.frame(spread_df)
row.names(spread_df)<-spread_df$isPaneth
spread_df$isPaneth<-NULL
print(spread_df)
chisq <- chisq.test(t(spread_df))
#dato_cet.summary_paneth$isPaneth<-as.factor(dato_cet.summary_paneth)
colori_ciclo<-c("#0173b4", "#93c2e9", "#b0d395")
ciclo<-ggplot(dato.summary_paneth,aes(x=" ",y=percent, fill=V2)) +
  geom_bar(width = 1, stat = "identity")+scale_fill_manual(values=colori_ciclo) +#geom_text(aes(y=mucca,label=ifelse(percent >= 0.03, paste0(sprintf("%.0f", percent*100),"%"),"")),colour="black")+
  coord_polar("y", start=0) +
  facet_grid(.~ isPaneth) +theme_void()+ggtitle(paste('chi_squared: pvalue:',round(chisq$p.value,15)))+theme(plot.title = element_text(hjust = 0.5))


design<-"
  1
  3
"
l<-ciclo+paneth+ plot_layout(design = design)

pdf(plot_out)
print(l)
graphics.off()


main <- function(log_f) {
  con <- file(log_f, open = "wt")
  sink(con)                   # cattura stdout (print, cat)
  sink(con, type = "message") # cattura stderr (message, warning)
  on.exit({
    flush(con)                # assicura che tutto sia scritto
    sink(type = "message")
    while (sink.number() > 0) sink()
    close(con)
  }, add = TRUE)

  cat("cellcycle percentage\n")
  print(dato.summary_paneth)

  cat("paneth percentage\n")
  print(dato.summary)
}


main(log_out)





