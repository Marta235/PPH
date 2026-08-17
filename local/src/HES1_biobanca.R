library(ggplot2)
library(ggrepel)

lmo_fpkm<-snakemake@input[['lmo_fpkm']]
lmx_fpkm<-snakemake@input[['lmx_fpkm']]
lmo_out<-snakemake@output[['out_lmo']]
lmx_out<-snakemake@output[['out_lmx']]
lmo_fpkm<-read.table(lmo_fpkm, header=T, sep=',')
lmx_fpkm<-read.table(lmx_fpkm, header=T, sep=',')
log_f <- snakemake@log[['log']]

  library(ggplot2)
  plot_nt_cet <- function( pd){
  nice_breaks_0 <- function(y, k = 5) {
    y_max <- max(y, na.rm = TRUE)
    
    if (!is.finite(y_max) || y_max <= 0) {
      return(list(
        breaks = c(0, 1),
        top = 1
      ))
    }
    
    raw_step <- y_max / k
    mag <- 10^floor(log10(raw_step))
    base <- raw_step / mag
    
    nice_base <- c(1, 2, 2.5, 5, 10)
    step <- nice_base[which(nice_base >= base)[1]] * mag
    
    top <- ceiling(y_max / step) * step
    breaks <- seq(0, top, by = step)
    
    list(
      breaks = breaks,
      top = top
    )
  }
  nice_axis_0 <- function(y, n = 5){
    y_max <- max(y, na.rm = TRUE)
    if (!is.finite(y_max) || y_max <= 0)
      return(list(breaks = c(0, 1), top = 1))
    
    raw_step <- y_max / n
    mag <- 10^floor(log10(raw_step))
    base <- raw_step / mag
    
    nice_base <- c(1, 2, 2.5, 5, 10)
    step <- nice_base[which(nice_base >= base)[1]] * mag
    
    top <- ceiling(y_max / step) * step
    breaks <- seq(0, top, by = step)
    
    list(breaks = breaks, top = top)
  }
  
  br <- nice_breaks_0(pd$FPKM, k = 5)  # scegli k fisso per tutti i plot
  
  p<-ggplot(pd, aes(x = Condition, y = FPKM)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0,size = 0.3) + theme_bw() +
    theme(
      text = element_text( size = 8),
      
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black")
    ) +
    scale_y_continuous(
      limits = c(0, br$top),
      breaks = br$breaks,
      expand = c(0, 0)
    )
  
  sink(log_f, append=T)
  print(nrow(pd))
  print(head(pd))
  sink()   
  print(t.test(formula=as.formula('FPKM~Condition'), data=pd))
  return(p)}

  ggsave(lmo_out, plot = plot_nt_cet(lmo_fpkm), width = 55,units = 'mm',height = 55)
  ggsave(lmx_out, plot = plot_nt_cet(lmx_fpkm), width = 55,units = 'mm',height = 55)
