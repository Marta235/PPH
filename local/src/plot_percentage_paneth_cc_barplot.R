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
cellpress_theme <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      # niente griglie
      panel.grid = element_blank(),

      # linee assi sottili
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      axis.ticks.length = unit(1.2, "mm"),

      # testo
      axis.text = element_text(size = base_size),
      axis.title = element_text(size = base_size),
      plot.title = element_text(size = base_size, hjust = 0.5, face = "plain"),

      # legenda (poi la sistemi in AI, ma intanto pulita)
      legend.title = element_blank(),
      legend.text = element_text(size = base_size),

      # margini compatti per figure piccole
      plot.margin = margin(1.5, 1.5, 1.5, 1.5, unit = "mm")
    )
}


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
colori_paneth <- c("Others"="#76069A", "Paneth"="#099963")

df_paneth <- dato.summary %>%
  mutate(
    value = ifelse(isPaneth == "Others", -percent, percent),
    isPaneth = factor(isPaneth, levels = c("Others", "Paneth"))
  )

colori_paneth <- c("Others" = "#76069A", "Paneth" = "#099963")

paneth<- ggplot(df_paneth, aes(x = isPaneth, y = value, fill = isPaneth)) +
  geom_col(width = 0.9) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(
    aes(
      label = paste0(round(percent * 100, 1), "%"),
      hjust = ifelse(value > 0, -0.15, 1.15)
    ),
    size = 2.6
  ) +coord_flip(clip = "off") +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.5),
    labels = function(x) paste0(abs(x) * 100, "%"),
    expand = expansion(mult = c(0.02, 0.10)) # spazio per testo fuori barra
  ) + 
  scale_x_discrete(
    expand = expansion(mult = c(0.2, 0.2))  # 👈 RIDUCE spazio verticale
  ) +
  scale_fill_manual(values = colori_paneth) +
  labs(x = NULL, y = NULL) +
  cellpress_theme(base_size = 6) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),legend.position = c(0.98, 0.98),
  legend.justification = c(1, 1),
  legend.background = element_blank(),
  legend.key = element_blank()
  )
paneth <- paneth+
  scale_x_discrete(expand = expansion(add = 0))+theme(
  legend.position = c(0.98, 0.98),
  legend.justification = c(1, 1),

  legend.direction = "vertical",
  legend.box = "vertical",

  legend.title = element_blank(),

  legend.text = element_text(size = 6),
  legend.key.size = unit(3, "mm"),

  legend.spacing.y = unit(1, "mm"),
  legend.spacing.x = unit(1, "mm"),

  legend.key = element_rect(fill = NA, colour = NA),
  legend.background = element_blank()
)



#cet
cont <- dato.summary_paneth[, c("isPaneth","V2","total_count")]
tab <- xtabs(total_count ~ isPaneth + V2, data = cont)   # 2x3
chisq <- chisq.test(tab)

# dati per plot (percentuale per gruppo)
df_plot <- cont %>%
  group_by(isPaneth) %>%
  mutate(percent = total_count / sum(total_count)) %>%
  ungroup() %>%
  mutate(
    value = ifelse(isPaneth == "Others", -percent, percent),   # diverging
    V2 = factor(V2, levels = c("G1","S","G2M")),               # ordine (cambialo se vuoi)
    isPaneth = factor(isPaneth, levels = c("Others","Paneth"))
  )

colori_ciclo <- c("G1"="#0173b4", "G2M"="#93c2e9", "S"="#b0d395")

ciclo <- ggplot(df_plot, aes(x = V2, y = value, fill = V2)) +
  geom_col(width = 0.9) +  
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(
    aes(
      label = paste0(round(percent * 100, 1), "%"),
      hjust = ifelse(value > 0, -0.12, 1.12)
    ),
    size = 2.6,      # ~7 pt in ggplot-ish (dipende dal device), ma a 56mm rimane leggibile
  
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.5),
    labels = function(x) paste0(abs(x) * 100, "%"),
    expand = expansion(mult = c(0.02, 0.10)) # spazio per testo fuori barra
  ) +
  scale_fill_manual(values = colori_ciclo) +
  labs(
    x = NULL,
    y = "Percent of cells"
    # title = NULL  # meglio in caption/AI
  ) +
  cellpress_theme(base_size = 6) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
  legend.position = c(0.98, 0.98),
  legend.justification = c(1, 1),
  legend.background = element_blank(),
  legend.key = element_blank(),
    plot.margin = margin(1.5, 6, 1.5, 1.5, unit = "mm")
  )

ciclo <- ciclo +
  scale_x_discrete(expand = expansion(add = 0))+theme(
  legend.position = c(0.98, 0.98),
  legend.justification = c(1, 1),

  legend.direction = "vertical",
  legend.box = "vertical",

  legend.title = element_blank(),

  legend.text = element_text(size = 6),
  legend.key.size = unit(3, "mm"),

  legend.spacing.y = unit(1, "mm"),
  legend.spacing.x = unit(1, "mm"),

  legend.key = element_rect(fill = NA, colour = NA),
  legend.background = element_blank()
)


design<-"
  1
  3
"
l <- paneth  / ciclo + plot_layout(heights = c(2, 3))



pdf(plot_out, width = 55/25.4, height = 55/25.4, useDingbats = FALSE)
print(l)
dev.off()


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





