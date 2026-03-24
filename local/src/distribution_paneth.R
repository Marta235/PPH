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
meta_path<-snakemake@input$metagene
clu_path<-snakemake@input$cluster



output_plot<-snakemake@output$out[1]
#kmeans<- read.table(file = clust,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
meta<- read.table(file = meta_path,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
meta<-as.data.frame(meta)

cluster<-read.table(file = clu_path,row.names = 1,sep=",",header = TRUE,stringsAsFactors = FALSE)
cluster<-as.data.frame(cluster)

df_combined<-merge(meta,cluster,by='row.names')
print(head(df_combined))
rownames(df_combined)<-df_combined$Row.names
df_combined$Row.names<-NULL
#df_combined<-df_combined[df_combined$isPaneth!='filtered',]
df_combined$isPaneth[df_combined$isPaneth=='Paneth']<-'SPC'
df_combined$isPaneth[df_combined$isPaneth=='Others']<-'NSPC'
df_combined$isPaneth[df_combined$isPaneth=='filtered']<-'NSPC'



x_min <- 0
x_max <- ceiling(max(df_combined$x))
breaks_x <- pretty(c(x_min, x_max), n = 5) 
x_max<-max(breaks_x)

df_combined$isPaneth <- factor(df_combined$isPaneth, levels = c("SPC", "NSPC"), labels = c("SPC", "NSPC"))

TOTCELL<-nrow(df_combined)
a <- ggplot(df_combined, aes(x = x, color = isPaneth)) +
  geom_density(aes(y = after_stat(count/TOTCELL)), position = "identity", bw = 0.05, size = 1)+
  scale_color_manual(name = "Cluster", 
                     values = c("NSPC" = "#76069A", "SPC" = "#099963")) +
  scale_x_continuous(expand = c(0, 0), limits = c(x_min, x_max), breaks = breaks_x) +  # Specifica i tick manualmente
  ggtitle("Metagene Distribution") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),  # Sfondo bianco
        panel.grid.major = element_blank(),                          # Rimuovi griglie maggiori
        panel.grid.minor = element_blank(),                          # Rimuovi griglie minori
        axis.line = element_line(color = "black"),                   # Colore nero per gli assi
        axis.ticks = element_line(color = "black"),                  # Tick marks neri
        axis.ticks.length = unit(0.2, "cm"),                         # Lunghezza dei tick
        axis.title.x = element_text(size = 12),                      # Etichetta asse X
        axis.title.y = element_text(size = 12))                      # Etichetta asse Y
  guides(color = guide_legend(override.aes = list(linetype = 1, size = 1, shape = NA, fill = NA)))


ggp <- ggplot_build(a)
data <- ggp$data[[1]]

data$y_orig <- data$y

maxy <- max(data$y_orig)
print('=============================================')
print(maxy)

breaks_y <- pretty(c(0, maxy), n = 4) 
maxy<-max(breaks_y)

data$coloripaneth <- ifelse(data$colour == "#76069A", 'NSPC', 'SPC')

TOTCELL<-nrow(df_combined)
a<-ggplot(df_combined, aes(x = x, color = isPaneth)) +
  geom_density(aes(y = after_stat(count/TOTCELL)), position = "identity", bw = 0.05,linewidth = 0.8,
    key_glyph = draw_key_path) +
  scale_color_manual(name = "Cluster", 
                     values = c("NSPC" = "#76069A", "SPC" = "#099963")) +
  scale_x_continuous(expand = c(0, 0), limits = c(x_min, x_max), breaks = breaks_x) + 
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),  # Sfondo bianco
        panel.grid.major = element_blank(),                          # Rimuovi griglie maggiori
        panel.grid.minor = element_blank(),                          # Rimuovi griglie minori
        axis.line = element_line(color = "black"),                   # Colore nero per gli assi
        axis.ticks = element_line(color = "black"),                  # Tick marks neri
        axis.ticks.length = unit(0.2, "cm"),                         # Lunghezza dei tick
        axis.title.x = element_text(size = 12),                      # Etichetta asse X
        axis.title.y = element_text(size = 12))  +                    # Etichetta asse Y
guides(color = guide_legend(override.aes = list(linetype = 1, size = 1, shape = NA, fill = NA)))+
        scale_y_continuous(expand = c(0, 0),breaks=breaks_y,limits=c(0, maxy)) +
  labs(y = "PDF* fraction of cells", x = "SPC Metagene")                      # Etichetta asse Y
  #guides(color = guide_legend(override.aes = list(linetype = 1, size = 1, shape = NA, fill = NA)))
a <- a + 
  theme_minimal() +
  theme(
  panel.background = element_rect(fill = "white", color = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.ticks.length = unit(0.2, "cm"),
  axis.title.x = element_text(size = 8, color = "black"),
  axis.title.y = element_text(size = 8, color = "black"),
  axis.text = element_text(color = "black"),
  legend.text = element_text(color = "black")
)+guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 1
      )
    )
  )

ggsave(output_plot, plot=a, width=100, height=100, units="mm")
 
#save.image(paste0(output_plot, '.Rdata'))
