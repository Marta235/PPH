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
#library(ggpubr)

#output_plot<-snakemake@output[[1]]
cet_data_path<-snakemake@input$geni_cetux
nt_data_path<-snakemake@input$geni_nt
cet_clu_path<-snakemake@input$cluster_cetux
nt_clu_path<-snakemake@input$cluster_nt



output_plot<-snakemake@output$out[1]
#kmeans<- read.table(file = clust,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)

cet<- read.table(file = cet_data_path,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
cet<-as.data.frame(cet)
cet$trattamento<-'Cetuximab'
cluster_cet<-read.table(file = cet_clu_path,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
cluster_cet<-as.data.frame(cluster_cet)

nt<- read.table(file = nt_data_path,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
nt<-as.data.frame(nt)
nt$trattamento<-'Untreated'
cluster_nt<-read.table(file = nt_clu_path,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
cluster_nt<-as.data.frame(cluster_nt) 

data<-rbind(cet,nt)
cluster<-rbind(cluster_cet,cluster_nt)

df_combined<-merge(data,cluster,by='row.names')
print(head(df_combined))
rownames(df_combined)<-df_combined$Row.names
df_combined$Row.names<-NULL
df_combined<-df_combined[df_combined$isPaneth!='filtered',]
df_combined<-df_combined[df_combined$isPaneth!='Paneth',]


# Combina i due dataframe


# Crea il grafico
x_min <- 0
x_max <- ceiling(max(df_combined$HES1))



breaks_x <- pretty(c(x_min, x_max), n = 5) 
x_max<-max(breaks_x)
print(breaks_x) 


a <- ggplot(df_combined, aes(x = HES1, fill = trattamento)) +
  geom_histogram( alpha=0.3, binwidth=0.05, position = 'identity')+
  
  #geom_density((aes(y=after_stat(scaled))),position = "identity", bw = 0.05, size = 0.8) +#
  scale_fill_manual(name = "Treatment", 
                     values = c("Cetuximab" = "red", "Untreated" = "black")) +
  scale_x_continuous(expand = c(0, 0), limits = c(x_min, x_max), breaks = breaks_x) +  # Specifica i tick manualmente
    # Mantiene l'asse Y gestito automaticamente
  labs(y = "Fraction of total cells(A.U.)", x = "HES1")+
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),  
        panel.grid.major = element_blank(),                          
        panel.grid.minor = element_blank(),                         
        axis.line = element_line(color = "black"),                   
        axis.ticks = element_line(color = "black"),                  
        axis.ticks.length = unit(0.2, "cm"),                         
        axis.title.x = element_text(size = 8),                      
        axis.title.y = element_text(size = 8),legend_position = 'none')                      
  guides(color = guide_legend(override.aes = list(linetype = 1, size = 1, shape = NA, fill = NA)))

density_data <- ggplot_build(a)$data[[1]]


maxy <- ceiling(max(density_data$y))# [[2]] if histogram is removed
breaks_y <- pretty(c(0, maxy), n = 5) 
maxy<-max(breaks_y)
print(maxy)
print(breaks_y)
a<-a+scale_y_continuous(expand = c(0, 0), limits = c(0, maxy),breaks=breaks_y) +
labs(y = "cells number", x = "HES1")+theme(,legend_position = 'none')
a <- a + theme(
  legend.position = c(1, 1),  # alto a destra
  legend.justification = c(1, 1)
)

print(head(density_data))

ggsave(output_plot, plot=a, width=100, height=100, units="mm")
