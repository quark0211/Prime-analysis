import os, sys
import numpy as np
import random
import anndata
import pandas as pd
import loompy as lp
import scanpy as sc

def nb_thresh(U,S,var_t = 1.5,u_min =0.02,s_min =0.02):
    '''
    Take in U,S matrices, and find genes that meet var/mean thresh
    U,S: cellxgene count matrices
    var_t: dispersion threshold
    u_min: minumum U average count
    s_min: minumum S average count
    '''
    var_threshold = var_t
    U_mean = U.mean(0)
    S_mean = S.mean(0)
    U_var = U.var(0)
    S_var = S.var(0)

    u_min = u_min
    s_min =  s_min


    fitted_idx = (U_mean > u_min) & (S_mean > s_min) \
    & (((U_var-U_mean)/(U_mean**2)) > var_threshold)\
    & (((S_var-S_mean)/(S_mean**2)) > var_threshold)\
    & (np.abs(np.log(S_mean/U_mean)) < 4)


    return fitted_idx


l = 'dataset/real_data/mouse_all.loom'
ds = lp.connect(l)
print(ds)
S = ds.layers['spliced'][:,:]
U = ds.layers['unspliced'][:,:]
barcode= ds.col_attrs['barcode']
true_labels = ds.ca['subclass_label']
g_names = ds.ra['gene_name']
print(barcode.shape)
print(S.shape)
ds.close()

#Filter for genes that pass thresholds
fitted_idx = nb_thresh(U.T,S.T,var_t = 1.5,u_min =0.02,s_min =0.02)
print('No. all genes that pass thresh: ', np.sum(fitted_idx))


#Get HVGs by standard methods
X=S.T  #May be better to use U.T for nuclear, snRNAseq, data
adata = anndata.AnnData(X=X)
adata.layers["counts"] = adata.X.copy()  # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
print(X.shape)


g_names_toUse = g_names[fitted_idx & adata.var.highly_variable]
print('No. of genes in 2k HVGs that pass filter: ',len(g_names_toUse))
print(g_names_toUse)

hvg_set = set(g_names_toUse)
indices = [i for i, gene in enumerate(g_names) if gene in hvg_set]

df_U = pd.DataFrame(U[indices,:])
df_S = pd.DataFrame(S[indices,:])
df_labels = pd.DataFrame(true_labels)
df_barcode = pd.DataFrame(barcode)

df_U.to_csv("mouse_hvg_2020_unspliced.csv")
df_S.to_csv("mouse_hvg_2020_spliced.csv")
df_labels.to_csv("mouse_subclass2020_labels.csv")