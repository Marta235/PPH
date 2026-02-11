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

#library(kSamples)
#test <- ks.test(cet, nt)
sink(snakemake@log[['log']])   
cat("KS test\n")
cat("======================\n\n")
#print(test)
sink()


x_min <- 0
#x_max <- ceiling(max(df_combined$x))
x_max<-12


breaks_x <- pretty(c(x_min, x_max), n = 5) 
x_max<-max(breaks_x)
print(breaks_x) 


orig <- ggplot(df_combined, aes(x = x, color = Treatment)) +
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

# we get the axis labels from the real/not rescaled densities plot
pdf_true_cet <- ggplot(df_combined[df_combined$Treatment=="Cetuximab",], aes(x = x, color = Treatment))+geom_density(bw = 0.05, size = 0.8)
pdf_true_nt <- ggplot(df_combined[df_combined$Treatment=="Not Treated",], aes(x = x, color = Treatment))+geom_density(bw = 0.05, size = 0.8)

density_data_nt_y <- ggplot_build(pdf_true_nt)$data[[1]]$y
density_data_cet_y <- ggplot_build(pdf_true_cet)$data[[1]]$y
density_scaleddata_y <- ggplot_build(orig)$data[[1]]$y

get_max_breaks <- function(my_d, digits=0) {
  maxy <- round(max(my_d), digits=digits)
  #breaks_y <- pretty(c(0, maxy), n = 5)  # pretty gives more than 5 for the real pdf axes
  breaks_y <- seq(0, maxy, length.out=6)
  maxy<-max(breaks_y)
  print(maxy)
  print(breaks_y)
  return(list(max=maxy, breaks=breaks_y))
}

pdf_left_nt_data <- get_max_breaks(density_data_nt_y, digits=1)
pdf_right_cet_data <- get_max_breaks(density_data_cet_y, digits=1)
scaled_pdf <- get_max_breaks(density_scaleddata_y)


a <- orig + scale_y_continuous(expand = c(0, 0), limits = c(0, scaled_pdf[['max']]),  breaks=scaled_pdf[['breaks']], labels=pdf_left_nt_data[['breaks']],
                               sec.axis=dup_axis(nam='PDF for Cetuximab Treated cells', labels=pdf_right_cet_data[['breaks']])) +
  labs(y = 'PDF for Not Treated cells', x = "Metagene")
a <- a + theme(
  legend.position = c(0.95, 0.95),  # alto a destra
  legend.justification = c(1, 1)
)

ggsave(output_plot, plot=a, width=100, height=100, units="mm")
# th<-6
# nt_dist<-density_data[density_data$colour=='black' ,'y']
# ctx_dist<-density_data[density_data$colour=='red','y']
# 
# ks.test(meta_nt[meta_nt$x>2,'x'], meta_cet[meta_cet$x>2,'x'])




#ggplot(ps, aes(x=x, y=y))+geom_point()



#save.image(paste0(output_plot, '.Rdata'))
