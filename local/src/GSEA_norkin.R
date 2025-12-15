### GSEA enrichment analysis

library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(ggplot2)

gene_res_f <- snakemake@input[["gene_res_freq"]]
GSEA_r <- snakemake@output[["GSEA_r"]]
GSEA_ridgeplot <- snakemake@output[["GSEA_ridgeplot"]]
#type <- snakemake@wildcards[["msign"]]
signature<-'Paneth'

gene_res_df <- read.table(gene_res_f, sep = "\t", header = TRUE,stringsAsFactors = FALSE)
score<-'Paneth_score'
###order
geneList <- gene_res_df[,score]
names(geneList) <- as.character(gene_res_df[,signature])
geneList <- sort(geneList, decreasing = TRUE)
print(head(geneList))



m_t2g <- read.table('/mnt/cold2/snaketree/prj/PPH/local/share/data/signature_paperPaneth/all_signature.csv', quote = "", sep = ",", header = FALSE,stringsAsFactors = FALSE )
#print(m_t2g)
print('==============================')
str(geneList)
head(names(geneList))

head(m_t2g)
str(m_t2g)

# quanta sovrapposizione c'è?
ol <- intersect(names(geneList), m_t2g[[2]])
length(ol)
head(ol)

em <- GSEA(geneList, TERM2GENE = m_t2g, pvalueCutoff = 1,nPerm=10000)

write.table(em@result, file = GSEA_r, quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)
pdf(GSEA_ridgeplot,width=12,height=12)
ridgeplot(em, showCategory = 15)
graphics.off()
ggsave(GSEA_ridgeplot)


