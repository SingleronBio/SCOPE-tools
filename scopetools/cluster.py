# -*- coding: utf-8 -*-
from scopetools.utils import getlogger, CommandWrapper
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np


sc.settings.set_figure_params(dpi=120)

if not os.path.exists(args['outdir']):
    os.mkdir(args['outdir'])
rg_method = args['rgenes_method']
hfile_name = '_' + args['input'].split("/")[-1]

# adata = sc.read(args['input'])
adata = sc.read_mtx(args['input'])
adata.var_names_make_unique()
os.chdir(args['outdir'])
sc.pl.highest_expr_genes(adata, n_top=20, save=hfile_name + ".pdf")
sc.pl.highest_expr_genes(adata, n_top=20, save=hfile_name + ".png")

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

mito_genes = adata.var_names.str.startswith("mt-")
if not mito_genes.all():
    mito_genes = adata.var_names.str.startswith("MT-")
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True, save=hfile_name + ".pdf")
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True, save=hfile_name + ".png")

# filter data
adata = adata[adata.obs.n_genes < args["max_genes"], :]
adata = adata[adata.obs.percent_mito < args["max_permito"], :]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, save=hfile_name + ".pdf")
sc.pl.highly_variable_genes(adata, save=hfile_name + ".png")

# Actually do the filtering
adata = adata[:, adata.var.highly_variable]
# 回归每个细胞总计数和线粒体基因表达百分比的影响。将数据放缩到方差为1。
sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
sc.pp.scale(adata, max_value=10)

# Principal component analysis
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, save=hfile_name + ".pdf")
sc.pl.pca_variance_ratio(adata, log=True, save=hfile_name + ".png")
sc.pp.neighbors(adata, n_neighbors=int(args['neighbors']), n_pcs=int(args['pc']))
# sc.tl.umap(adata)
# ngenes=len(adata.var)
if args['cluster_algo'] == "leiden":
    # sc.settings.verbosity=2
    sc.tl.leiden(adata)
    sc.tl.rank_genes_groups(adata, 'leiden', method=rg_method, n_genes=100000)
# sc.pl.rank_genes_groups(adata,n_genes=25,sharey=False,save=hfile_name+".png")
elif args['cluster_algo'] == "louvain":
    sc.tl.louvain(adata)
    sc.tl.rank_genes_groups(adata, 'louvain', method=rg_method, n_genes=100000)

if args['plot_method'] == "tsne":
    sc.tl.tsne(adata)
    sc.pl.tsne(adata, color=[args['cluster_algo']], save=hfile_name + ".pdf")
    sc.pl.tsne(adata, color=[args['cluster_algo']], save=hfile_name + ".png")
elif args['plot_method'] == "umap":
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=[args['cluster_algo']], save=hfile_name + ".pdf")
    sc.pl.umap(adata, color=[args['cluster_algo']], save=hfile_name + ".png")

result = adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names
s = pd.DataFrame(
    {str(group) + '_' + key: result[key][group]
     for group in groups for key in ['names', 'logfoldchanges', 'pvals', 'pvals_adj']})
for i in groups:
    file = "cluster" + i + "_diffgenes.csv"
    s.iloc[:, int(i) * 4:(int(i) + 1) * 4].to_csv(file, index=None, header=['names', 'logFC', 'pvals', 'pvals_adj'])

adata.obs.to_csv("./obs.csv", index=None, header=True)
adata.var.to_csv("./var.csv", index=None, header=True)
adata.write("./results.h5ad")


def cluster(ctx, input, outdir, filter_genome, max_genes, max_permito, neighbors, pc, rgenes_method, cluster_algo, plot_method):
    sc.settings.set_figure_params(dpi=120)
    pass
