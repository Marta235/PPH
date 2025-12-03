library(data.table)
library(mclust)
library(cluster)
library(reshape2)
library(ggplot2)
library(ramify)
library(psych)
library(pheatmap)
library(resample)
library(dplyr)
library(viridis)
library(igraph)
library(patchwork)

input      <- snakemake@input[['data']]
kmeans_out <- snakemake@output[['out']]
log_f      <- snakemake@log[['log']]

main <- function() {
  # apri sink per stdout e stderr (message/warning)
  con <- file(log_f, open = "wt")
  sink(con)                          # stdout
  sink(con, type = "message")        # stderr
  on.exit({
    sink(type = "message")
    while (sink.number() > 0) sink()
    close(con)
  }, add = TRUE)

  dato   <- read.table(file = input, row.names = 1, sep = ",", header = TRUE)
  dato_t <- transpose(dato)
  rownames(dato_t) <- colnames(dato)
  colnames(dato_t) <- rownames(dato)

  stringa_split <- function(stringa) strsplit(stringa, split = ":")[[1]][2]
  colnames(dato_t) <- sapply(colnames(dato_t), FUN = stringa_split)

  # Select Paneth cell markers
  cinque <- c("ATOH1","GFI1","DLL1","DEFA5","DEFA6")

  presenti <- intersect(cinque, colnames(dato_t))
  mancanti <- setdiff(cinque, presenti)

  if (length(mancanti) > 0) {
    message("Geni non trovati e rimossi: ", paste(mancanti, collapse = ", "))
  }
  if (length(presenti) == 0) {
    message("Nessuno dei geni richiesti è presente in dato_t. Creo file vuoto.")
    dir.create(dirname(kmeans_out), recursive = TRUE, showWarnings = FALSE)
    file.create(kmeans_out)
    return(invisible(NULL))   # esci pulito: on.exit chiuderà i sink
  }

  cinque_df <- dato_t[, presenti, drop = FALSE]


  set.seed(123)
  cl <- kmeans(cinque_df, centers = 2)

  meta_mu <- apply(cl$centers, 1, geometric.mean)
  ordine  <- order(unlist(meta_mu))
  ordinate <- c('Others','Paneth')

  final_ordinate <- c()
  for (i in seq(1,2)) final_ordinate[ordine[i]] <- ordinate[i]

  cluster_id <- cl$cluster
  isPaneth <- vapply(cluster_id, function(el) final_ordinate[el], character(1))
  posteriors <- data.frame(cluster_id = as.numeric(cluster_id), isPaneth = isPaneth)

  conteggio <- table(posteriors$isPaneth)
  print(conteggio)

  write.table(posteriors, file = kmeans_out, sep = ',', quote = FALSE, row.names = FALSE)
}
main <- function() {

  con <- file(log_f, open = "wt")
  sink(con)                          # stdout
  sink(con, type = "message")        # stderr
  on.exit({
    sink(type = "message")
    while (sink.number() > 0) sink()
    close(con)
  }, add = TRUE)

main()





