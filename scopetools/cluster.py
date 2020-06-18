# -*- coding: utf-8 -*-
import os

import pandas as pd
import scanpy as sc

from .report import Reporter
from .utils import getlogger

FILETYPE = ['png', 'pdf']
CLUSTER_ALGORITHM = ['leiden', 'louvain']
sc.settings.set_figure_params(dpi_save=150)
sc.settings.verbosity = 0

logger = getlogger(__name__)
logger.setLevel(10)


def cluster(ctx, matrix, outdir, sample, barcodes, genes, n_top, min_genes, min_cells, n_genes_by_counts, pct_counts_mt, exclude_highly_expressed, max_fraction, n_top_genes, max_value, n_neighbors, n_pcs):
    sample_outdir = outdir / sample / '06.cluster'
    sample_outdir.mkdir(parents=True, exist_ok=True)
    os.chdir(sample_outdir)

    logger.info('cluster start!')

    adata_mtx = sc.read_mtx(matrix).T
    obs = pd.read_csv(barcodes, index_col=0, header=None)
    var = pd.read_csv(genes, index_col=0, header=None)
    obs.index.set_names(None, inplace=True)
    var.index.set_names(None, inplace=True)
    adata_mtx.obs = obs
    adata_mtx.var = var

    result_file = sample_outdir / f'{sample}.h5ad'

    adata_mtx.var_names_make_unique()
    [sc.pl.highest_expr_genes(adata_mtx, n_top=n_top, save=f'_{sample}.{filetype}') for filetype in FILETYPE]

    sc.pp.filter_cells(adata_mtx, min_genes=min_genes)
    sc.pp.filter_genes(adata_mtx, min_cells=min_cells)

    adata_mtx.var['mt'] = adata_mtx.var_names.str.upper().str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata_mtx, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    [sc.pl.violin(adata_mtx, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save=f'_{sample}.{filetype}') for filetype in FILETYPE]
    [sc.pl.scatter(adata_mtx, x='total_counts', y='pct_counts_mt', save=f'_pct_mt_{sample}.{filetype}') for filetype in FILETYPE]
    [sc.pl.scatter(adata_mtx, x='total_counts', y='n_genes_by_counts', save=f'_n_genes_{sample}.{filetype}') for filetype in FILETYPE]

    adata_mtx = adata_mtx[adata_mtx.obs.n_genes_by_counts < n_genes_by_counts, :]
    adata_mtx = adata_mtx[adata_mtx.obs.pct_counts_mt < pct_counts_mt, :]
    sc.pp.normalize_total(adata_mtx, target_sum=1e6, exclude_highly_expressed=exclude_highly_expressed, max_fraction=max_fraction)
    sc.pp.log1p(adata_mtx)
    sc.pp.highly_variable_genes(adata_mtx, n_top_genes=n_top_genes)

    [sc.pl.highly_variable_genes(adata_mtx, save=f'_{sample}.{filetype}') for filetype in FILETYPE]

    adata_mtx.raw = adata_mtx
    adata_mtx = adata_mtx[:, adata_mtx.var.highly_variable]

    sc.pp.regress_out(adata_mtx, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata_mtx, max_value=max_value)
    sc.tl.pca(adata_mtx, svd_solver='arpack')

    sc.pp.neighbors(adata_mtx, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata_mtx)
    sc.tl.leiden(adata_mtx)
    sc.tl.louvain(adata_mtx)
    [sc.pl.umap(adata_mtx, color=algo, save=f'_{algo}_{sample}.{filetype}') for filetype in FILETYPE for algo in CLUSTER_ALGORITHM]

    adata_mtx.write(result_file, compression='gzip')

    logger.info('cluster done!')

    stat_info = {}
    img = {}
    pngs = (sample_outdir / 'figures').rglob('*.png')
    for png in pngs:
        img[png.name] = png.resolve()

    # report
    logger.info('generate report start!')
    Reporter(name='cluster', stat_json=stat_info, outdir=sample_outdir.parent, img=img)
    logger.info('generate report done!')
