library(ggplot2)
library(ggrepel)

set.seed(42)
## lista sensibili da PDX definitiva vol < -50% a 3 settimane
cetuxi <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/treats/last_march2024/last_cet_march2024.txt', sep="\t", header=TRUE, stringsAsFactors = F)
colnames(cetuxi)[2] <- 'Cetuxi_VolVar3WKS'
colnames(cetuxi)[3] <- 'Cetuxi_VolVar6WKS'
cetuxi[2] <- cetuxi[2]*100
cetuxi[3] <- cetuxi[3]*100

sens <- cetuxi[cetuxi$Cetuxi_VolVar3WKS <= -50,]

## lista 4wt no MET no HER
genetic_res <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/genetic_cet_res_update0625.tsv', sep="\t", header=FALSE, stringsAsFactors = F)
colnames(genetic_res) <- 'CASE'

## lista badboys
badboys <- read.table('/mnt/trcanmed/snaketree/prj/pdxopedia/local/share/data/badboys', sep="\t", header=FALSE, stringsAsFactors = F)
colnames(badboys) <- 'CASE'

universe <- setdiff(sens$CASE, genetic_res$CASE)
good_universe <- setdiff(universe, badboys$CASE)

length(universe)
length(good_universe)

# 48 sensibili, non ci sono mutati

## {lista PDO validati} per ora no
#pdo <- read.table('/mnt/trcanmed/snaketree/prj/biobanca/local/share/data/whoiswho_validation_xen_revision_derivation.tsv', sep="\t", header=TRUE, stringsAsFactors = F)
#colnames(pdo)[2] <- 'PDO Validation'
#colnames(pdo)[3] <- 'PDO Derivation'

# metadata rnaseq
rna <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/samples_data', sep="\t", stringsAsFactors = F, header=T)

treatments <- unique(rna$type[grepl('cetux', rna$type)])

rna_keep <- rna[rna$type %in% treatments, ]
rna_keep$CASE <- substr(rna_keep$id, 0, 7)

rna_keep_2 <- rna_keep[rna_keep$CASE %in% good_universe,]

nrow(rna_keep)
nrow(rna_keep_2)
length(unique(rna_keep_2$CASE))

length(unique(rna_keep_2$CASE[grepl('LMX', rna_keep_2$type)]))
length(unique(rna_keep_2$CASE[grepl('LMO', rna_keep_2$type)]))
# 33 con dato di expr, 33 pdx 12 pdo (ok in linea con DEG storici)
lmo_sel <- rna_keep_2[grepl('LMO', rna_keep_2$type),]
lmx_sel <- rna_keep_2[grepl('LMX', rna_keep_2$type),]

fpkm <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/fpkm_H.tsv.gz'), sep="\t", stringsAsFactors = F, header=T)

get_ave <- function(model, data) {
  myd <- data[, grepl(model, colnames(data)), drop=F]
  rowMeans(myd)
}

# we do not have .1 issues in cetuxi experiments
plot_nt_cet <- function(cases, expr, gene='HES1') {
  nt <- cases[grepl('NT', cases$type),]
  cet <- cases[!grepl('NT', cases$type),]
  sel_expr_nt <- expr[, colnames(expr) %in% nt$id, ]
  sel_expr_cet <- expr[, colnames(expr) %in% cet$id, ]
  cetm <- unique(cet$CASE)
  ntm <- unique(nt$CASE)
  print(setdiff(cetm, ntm))
  print(setdiff(ntm, cet))
  ave_nt <- sapply(ntm, get_ave, sel_expr_nt)
  ave_cet <- sapply(cetm, get_ave, sel_expr_cet)
  pd <- data.frame(FPKM=c(ave_nt[gene, ], ave_cet[gene,]), Condition=c(rep('NT', ncol(ave_nt)), rep('cetuximab', ncol(ave_cet))))
                   #model=c(colnames(ave_nt), colnames(ave_cet)))
  pd$Condition <- factor(pd$Condition, levels=c('NT', 'cetuximab'))
  print(ggplot(data=pd, aes(x=Condition, y=FPKM))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0)+theme_bw(base_size = 20))
  return(t.test(formula=as.formula('FPKM~Condition'), data=pd))
}

plot_nt_cet(lmx_sel, fpkm)
plot_nt_cet(lmo_sel, fpkm)

# provare con FPKM dei due sottogruppi
fpkm_o <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDO_72h_S/fpkm_H.tsv.gz'), sep="\t", stringsAsFactors = F, header=T)
fpkm_x <- read.table(gzfile('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5starOK_cetuxi_treat_PDX_S/fpkm_H.tsv.gz'), sep="\t", stringsAsFactors = F, header=T)


plot_nt_cet(lmx_sel, fpkm_x)
plot_nt_cet(lmo_sel, fpkm_o)
# provare 6W/72h separate per cet
lmx_72h <- lmx_sel[grepl('72h', lmx_sel$type),]
lmx_6w <- lmx_sel[!grepl('72h', lmx_sel$type),]

plot_nt_cet(lmx_72h, fpkm)
plot_nt_cet(lmx_6w, fpkm)

## farne uno con linee paired ? ma ci sono buchi, nope

###
sixwid <- nt$id[grepl('6w', nt$type)]
sixwidc <- cet$id[grepl('6w', cet$type)]
long <- unique(substr(sixwid, 0, 7), substr(sixwidc, 0, 7))
pd$is_sixw <- ifelse(c(colnames(ave_nt), colnames(ave_cet)) %in% long, 'yes', 'no')
table(pd$is_sixw)
print(ggplot(data=pd, aes(x=Condition, y=FPKM))+geom_boxplot(outlier.shape=NA)+geom_jitter(aes(color=is_sixw), height=0)+theme_bw(base_size = 20))
pd[pd$is_sixw=='yes',]



## logFC


ave_log_fc <-function(model, wexpr, gene, ss, pc=0.1) { # 0.0001 or 1, 1 was used for pis
  myss <- ss[grepl(model, ss$id),]
  mywexpr <- t(wexpr[,colnames(wexpr) %in% myss$id])
  m <- merge(myss, mywexpr, by.x='id', by.y="row.names")
  # we want at least 1 NT and 1 cet
  if (!(any(grepl('NT', m$type)) && any(!grepl('NT', m$type)))) {
    return(NA)
  }
  rep1 <- m[grepl('.1', m$type, fixed=T),]
  rep0 <- m[!grepl('.1', m$type, fixed=T),]
  lfcs <- c()
  if (nrow(rep1) == 2) {
    nt <- rep1[grepl('NT', rep1$type), gene]
    cet <- rep1[!grepl('NT', rep1$type), gene]
    lfcs <- c(lfcs, log2((cet+pc)/(nt+pc)))
  } 
  if (nrow(rep0) == 2) {
    nt <- rep0[grepl('NT', rep0$type), gene]
    cet <- rep0[!grepl('NT', rep0$type), gene]
    lfcs <- c(lfcs, log2((cet+pc)/(nt+pc)))
  }
  if (length(lfcs)>=1) {
    return(mean(lfcs))
  } else {
    return(NA) # we have unpaired rep
  }
}


two_genes <- function(cases, expr, geney, genex) {
  wexpr <- expr[rownames(expr) %in% c(genex, geney),]
  models <- unique(substr(cases$id, 0,7))
  
  lfcs_y <- sapply(models, ave_log_fc, wexpr, geney, cases)
  lfcs_x <- sapply(models, ave_log_fc, wexpr, genex, cases)
  return(data.frame(x=lfcs_x, y=lfcs_y, models=models))
}

xeno <- two_genes(lmx_sel, fpkm, 'HES1', 'DLL1')

xeno <- xeno[!is.na(xeno$x) & !is.na(xeno$y),]

ggplot(data=xeno, aes(x=x, y=y))+geom_point()+geom_smooth(method='lm')+theme_bw(base_size = 20)+xlab('DLL1')+ylab('HES1')+
  geom_text_repel(aes(label = models))


cor.test(xeno$x,xeno$y)


lmo <- two_genes(lmo_sel, fpkm, 'HES1', 'DLL1')

lmo <- lmo[!is.na(lmo$x) & !is.na(lmo$y),]

ggplot(data=lmo, aes(x=x, y=y))+geom_point()+geom_smooth(method='lm')+theme_bw(base_size = 20)+xlab('DLL1')+ylab('HES1')+
  geom_text_repel(aes(label = models))

cor.test(lmo$x,lmo$y)


###
for (g in c('ATOH1', 'GFI1', 'DEFA5', 'DEFA6')) {

xeno <- two_genes(lmx_sel, fpkm, 'HES1', g)

xeno <- xeno[!is.na(xeno$x) & !is.na(xeno$y),]

print(ggplot(data=xeno, aes(x=x, y=y))+geom_point()+geom_smooth(method='lm')+theme_bw(base_size = 20)+xlab(g)+ylab('HES1')+
  geom_text_repel(aes(label = models)))

print(cor.test(xeno$x,xeno$y))


lmo <- two_genes(lmo_sel, fpkm, 'HES1', g)

lmo <- lmo[!is.na(lmo$x) & !is.na(lmo$y),]

ggplot(data=lmo, aes(x=x, y=y))+geom_point()+geom_smooth(method='lm')+theme_bw(base_size = 20)+xlab(g)+ylab('HES1')+
  geom_text_repel(aes(label = models))

print(cor.test(lmo$x,lmo$y))

}