library(circlize)
library(ComplexHeatmap)
# xeno cet vs NT (S defined as OR in vivo, no genetic resistances here):
#gsea_path <- '/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDX_S/GSEA_results_C2_treat_cutoff0.05-cetuxi.vs.NT.tsv'
# pdo cet vs NT (S defined as OR in vivo, no genetic resistances here):
#gsea_path_cloni <- '/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/GSEA_results_C2_treat_cutoff0.05-cetuxi.vs.NT.tsv'
# xeno 4wt/noampl OR vs PD:
gsea_path <- '/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/cetuxi_resp_4wt_PDX/GSEA_results_C2_cetuxi_cutoff0.05-res.vs.sens.tsv'
# pdo 4wt/noampl OR vs PD (like before response defined in vivo): 
gsea_path_cloni <- '/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/cetuxi_resp_4wt_PDO/GSEA_results_C2_cetuxi_cutoff0.05-res.vs.sens.tsv'

# sfrutto codice di Sofia per mostrare Xeno e Pdo insieme, Xeno prendono il posto dei cas9 e Pdo dei cloni
data_gsea<-read.table(gsea_path,header = TRUE,sep='\t')
notch_gsea<-data_gsea[grep('NOTCH',data_gsea$ID),]
#notch_gsea<-notch_gsea[notch_gsea$enrichmentScore>0,]
#con_hes<-notch_gsea[grep('HES1',notch_gsea$core_enrichment),]

data_gsea_cloni<-read.table(gsea_path_cloni,header = TRUE,sep='\t')
#notch_gsea_cloni<-notch_gsea_cloni[notch_gsea_cloni$enrichmentScore>0,]
#con_hes_cloni<-notch_gsea_cloni[grep('HES1',notch_gsea_cloni$core_enrichment),]

per_heatmap<-notch_gsea[,c('ID','NES','pvalue','p.adjust')]
colnames(per_heatmap)<-paste0(colnames(per_heatmap),'_PDX')
hes_cloni<-data_gsea_cloni[rownames(per_heatmap),]
hes_cloni$NES_PDO=hes_cloni$NES
hes_cloni$pval_PDO<-hes_cloni$pvalue
hes_cloni$p.adjust_PDO<-hes_cloni$p.adjust
hes_cloni<-hes_cloni[,c('ID','NES_PDO','pval_PDO','p.adjust_PDO')]
per_heatmap<-merge(hes_cloni,per_heatmap,by="row.names")

rownames(per_heatmap)<-per_heatmap$Row.names
per_heatmap$Row.names<-NULL
per_heatmap$ID_PDX<-NULL
per_heatmap$ID_PDO<-NULL

nes_cols <- c("NES_PDX", "NES_PDO")
pval_cols <- setdiff(colnames(per_heatmap), nes_cols)

mat_nes <- per_heatmap[, nes_cols]
mat_pval <- per_heatmap[, pval_cols]

mat_nes<-as.matrix(mat_nes)
limite <- max((mat_nes), na.rm = TRUE)
#colorRampPalette(c("blue","white","red"))(100)
#limite<-1
colori<-col_fun <- colorRamp2(
  c(-limite, 0, limite),
  c("#4575b4", "white", "#d73027")  
)

ht_nes <- Heatmap(mat_nes,
                  name = "NES",
                  col = colori,
                  cluster_columns = FALSE,
                  cluster_rows = TRUE,
                  show_row_dend = FALSE ,
                  show_row_names = TRUE,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 8),
                  row_names_max_width = unit(12, "cm"),
                  heatmap_legend_param = list(
                    at = c(-limite, 0, limite),
                    labels = c(paste0("-", round(limite, 2)), "0", paste0("+", round(limite, 2))))
)
#print(ht_nes)

mat_pval<-as.matrix(mat_pval)
pval_th<-0.2 # ricapire un attimo questo...rifacciamo padjust solo su signature NOTCH o ...?
mat_pval_cat <- ifelse(mat_pval < pval_th, "significativo", "non significativo")
mat_pval_cat[is.na(mat_pval)] <- "NA" 

mat_pval_cat <- as.matrix(mat_pval_cat)
ordered_columns <- c(
  grep("_PDX$", colnames(mat_pval_cat), value = TRUE),
  grep("_PDO$", colnames(mat_pval_cat), value = TRUE)
)


mat_pval_cat <- mat_pval_cat[, ordered_columns]

colori_pval <- c("significativo" = "goldenrod", "non significativo" = "lightgrey", "NA" = "white")


ht_pval <- Heatmap(mat_pval_cat,
                   name = "p-value",
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   show_row_names = FALSE,
                   col = colori_pval,
                   rect_gp = gpar(col = "white", lwd = 1),
                   heatmap_legend_param = list(
                     at = c("significativo", "non significativo", "NA"),
                     labels = c(paste0("<", pval_th), paste0("≥ ", pval_th), "NA")
                   ))

draw(ht_nes + ht_pval, heatmap_legend_side = "right")