library(data.table)
library(ggplot2)
library(uwot)


cetux_input<-snakemake@input[['cetux']]
nt_input<-snakemake@input[['nt']]
log_out<-snakemake@output[['out']]
out_folder<-'projected_umap'
print(snakemake@wildcards[['sample_cet']])
cetux_out<-paste0(paste(out_folder,snakemake@wildcards[['sample_cet']],sep='/'),'_umap.csv')
nt_out<-paste0(paste(out_folder,snakemake@wildcards[['sample_nt']],sep='/'),'_umap.csv')

print(cetux_out)
print(nt_out)

min_dist<-0.3
# n_neighbors<-50

cetux<- read.table(file = cetux_input,row.names = 1,sep=",",header = TRUE)
cetux_t<- transpose(cetux)
rownames(cetux_t) <- colnames(cetux)
colnames(cetux_t)<-rownames(cetux)

nt<- read.table(file = nt_input,row.names = 1,sep=",",header = TRUE)
nt_t<- transpose(nt)
rownames(nt_t) <- colnames(nt)
colnames(nt_t)<-rownames(nt)

min_dist <- 0.8
n_neighbors <- 70 

set.seed(42)

umap_model <- umap(
  cetux_t,
  n_neighbors = n_neighbors,
  min_dist = min_dist,
  ret_model = TRUE
)


umap_cetux <- as.data.frame(umap_model$embedding)
colnames(umap_cetux) <- c("x", "y")
write.csv(umap_cetux, cetux_out, row.names = TRUE)


umap_nt <- umap_transform(nt_t, umap_model)
umap_nt <- as.data.frame(umap_nt)
colnames(umap_nt) <- c("x", "y")
write.csv(umap_nt, nt_out, row.names = TRUE)

main <- function(log_f) {

  con <- file(log_f, open = "wt")
  sink(con)                          # stdout
  sink(con, type = "message")        # stderr
  on.exit({
    sink(type = "message")
    while (sink.number() > 0) sink()
    close(con)
  }, add = TRUE)
    print("UMAP projection completed successfully.")
  } 



main(log_out)