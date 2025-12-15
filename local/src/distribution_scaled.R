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
library(diptest)
#library(ggpubr)

output_plot<-snakemake@output[[1]]
cet<-snakemake@input$c[1]
#cet<-'/mnt/cold1/snaketree/prj/scRNA/dataset/KMeans/kmeans/metageni/CRC0327_cetux_2_metagene.csv'
nt<-snakemake@input$c[2]
#nt<-'/mnt/cold1/snaketree/prj/scRNA/dataset/KMeans/kmeans/metageni/CRC0327_NT_2_metagene.csv'
#output_plot<-snakemake@output$out[1]

meta_cet<- read.table(file = cet,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
meta_cet<-as.data.frame(meta_cet)
meta_nt<-read.table(file = nt,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
meta_nt<-as.data.frame(meta_nt)
print(length(meta_cet$x))
print(length(meta_nt$x))
meta_nt$Treatment <- "Not Treated"
meta_cet$Treatment<- "Cetuximab"



df_combined <- rbind(meta_cet, meta_nt)
cet<-meta_cet$x
nt<-meta_nt$x

library(kSamples)
test <- ks.test(cet, nt)
sink(snakemake@log[['log']])   
cat("KS test\n")
cat("======================\n\n")
print(test)
sink()


x_min <- 0
x_max <- ceiling(max(df_combined$x))


breaks_x <- pretty(c(x_min, x_max), n = 5) 
x_max<-max(breaks_x)
print(breaks_x) 


a <- ggplot(df_combined, aes(x = x, color = Treatment)) +
  #geom_histogram(aes(fill='white'), alpha=0.3, binwidth=0.05, position = 'identity')+
  
  geom_density((aes(y=after_stat(scaled))),position = "identity", bw = 0.05, size = 0.8) +#
  scale_color_manual(name = "Treatment", 
                     values = c("Cetuximab" = "red", "Not Treated" = "black")) +
  scale_x_continuous(expand = c(0, 0), limits = c(x_min, x_max), breaks = breaks_x) +  # Specifica i tick manualmente
    # Mantiene l'asse Y gestito automaticamente
  labs(y = "Fraction of total cells(A.U.)", x = "Metagene")+
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
labs(y = "Relative cells number (A.U.)", x = "Metagene")+theme(,legend_position = 'none')
a <- a + theme(
  legend.position = c(0.95, 0.95),  # alto a destra
  legend.justification = c(1, 1)
)

print(head(density_data))

ggsave(output_plot, plot=a, width=100, height=100, units="mm")
# th<-6
# nt_dist<-density_data[density_data$colour=='black' ,'y']
# ctx_dist<-density_data[density_data$colour=='red','y']
# 
# ks.test(meta_nt[meta_nt$x>2,'x'], meta_cet[meta_cet$x>2,'x'])




#ggplot(ps, aes(x=x, y=y))+geom_point()



