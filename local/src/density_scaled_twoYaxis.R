library(ggplot2)

load('/mnt/cold2/snaketree/prj/PPH/dataset/Kmeans/distribution_scaled/CRC0322_cetux_1vsCRC0322_NT_1_3000_dist.pdf.Rdata')

 th <-  theme_minimal() +
   theme(panel.background = element_rect(fill = "white", color = NA),  
         panel.grid.major = element_blank(),                          
         panel.grid.minor = element_blank(),                         
         axis.line = element_line(color = "black"),                   
         axis.ticks = element_line(color = "black"),                  
         axis.ticks.length = unit(0.2, "cm"),                         
         axis.title.x = element_text(size = 8),                      
         axis.title.y = element_text(size = 8),legend_position = 'none') 
# 
# ggplot(df_combined, aes(x = x, fill = Treatment))+geom_histogram(alpha=0.3, binwidth=0.05, position = 'identity')+
#   scale_fill_manual(name = "Treatment", 
#                     values = c("Cetuximab" = "red", "Not Treated" = "black")) +
#   th
#          
# df_combined$x2 <- df_combined$x*0.1
# 
# pdf_true <- ggplot(df_combined, aes(x = x, color = Treatment))+geom_density(bw = 0.05, size = 0.8)+
#   scale_color_manual(name = "Treatment", 
#                     values = c("Cetuximab" = "red", "Not Treated" = "black")) +
#   th+            
#   guides(color = guide_legend(override.aes = list(linetype = 1, size = 1, shape = NA, fill = NA)))


pdf_true_cet <- ggplot(df_combined[df_combined$Treatment=="Cetuximab",], aes(x = x, color = Treatment))+geom_density(bw = 0.05, size = 0.8)+
  th+            
  guides(color = guide_legend(override.aes = list(linetype = 1, size = 1, shape = NA, fill = NA)))

pdf_true_nt <- ggplot(df_combined[df_combined$Treatment=="Not Treated",], aes(x = x, color = Treatment))+geom_density(bw = 0.05, size = 0.8)+
  th+            
  guides(color = guide_legend(override.aes = list(linetype = 1, size = 1, shape = NA, fill = NA)))

#pdf_true_nt <- ggplot(df_combined, aes(x = x2, color = Treatment))+geom_density(bw = 0.05, size = 0.8)+
#  scale_color_manual(name = "Treatment", 
#                     values = c("Cetuximab" = "red", "Not Treated" = "black")) +
#  th+            
#  guides(color = guide_legend(override.aes = list(linetype = 1, size = 1, shape = NA, fill = NA)))



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

###
load('/mnt/cold2/snaketree/prj/PPH/dataset/Kmeans/distribution_paneth/CRC0322_cetux_1_paneth_dist.pdf.Rdata')
ggplot(df_combined, aes(x = x, color = isPaneth)) +
  geom_density(position = "identity", bw = 0.05, size = 1) +
  scale_color_manual(name = "Cluster", 
                     values = c("Others" = "#76069A", "PCL cells" = "#099963")) +
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

ggplot(df_combined, aes(x = x, color = isPaneth)) +
  geom_density(aes(y = after_stat(count)), position = "identity", bw = 0.05, size = 1) +
  scale_color_manual(name = "Cluster", 
                     values = c("Others" = "#76069A", "PCL cells" = "#099963")) +
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


ggplot(df_combined, aes(x = x, fill = isPaneth))+geom_histogram(alpha=0.3, binwidth=0.05, position = 'identity')+
  scale_fill_manual(name = "Cluster", 
                     values = c("Others" = "#76069A", "PCL cells" = "#099963"))+th