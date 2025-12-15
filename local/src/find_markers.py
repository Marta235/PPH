import pandas as pd
import scanpy as sc
import seaborn as sns
import os

log_2data=snakemake.input.data
clus_data=snakemake.input.cluster
#leggo input e mergio
data=pd.read_csv(log_2data,index_col=0,header=0).T
gtf_dict={}
for col in data.columns:
    gtf_dict[col]=col.split(':')[1]
data.rename(columns=gtf_dict,inplace=True)

cluster_data=pd.read_csv(clus_data)

cluster_data=cluster_data.loc[:,'isPaneth']

data=pd.merge(data,cluster_data,left_index=True,right_index=True)

data = data.drop(data[data.isPaneth=='filtered'].index)


# ciclo su tutti i file nella cartella salvo i path che contengono nome 'gene_median', 'data'

sc_ = sc.AnnData(data.drop(columns='isPaneth'))
sc_.obs['isPaneth'] = pd.Categorical(data.loc[:, 'isPaneth'])
sc.tl.rank_genes_groups(sc_, 'isPaneth',method='wilcoxon',kwds=['only_positive'] )

nome_png='_'+snakemake.wildcards.sample+'.pdf'
a=sc.pl.rank_genes_groups(sc_,fontsize=7,save=nome_png,show=False, n_genes=40)
a=pd.DataFrame(sc_.uns['rank_genes_groups']['names'])
score=pd.DataFrame(sc_.uns['rank_genes_groups']['scores'])
score=score.rename(columns={'Paneth':'Paneth_score','Others':'Others_score'})
print(score.head())
res=pd.concat([a, score], axis=1)
print(res.head())
paneth=res.loc[:,['Paneth','Paneth_score']]

non_paneth=res.loc[:,['Others','Others_score']]
paneth.to_csv(snakemake.output.marker_paneth,sep='\t')
#non_paneth.to_csv(snakemake.output.marker_nonpaneth,sep='\t')

#pd.DataFrame(sc_.uns['rank_genes_groups']['pvals_adj']).to_csv(nome_csv+'pvals_adj.csv')


