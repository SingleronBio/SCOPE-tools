# -*- coding: utf-8 -*-
import os
from pathlib import Path

import pandas as pd
import scanpy as sc
import json
from scopetools.report import Reporter
from scopetools.utils import getlogger

FILETYPE = ['png', 'pdf']
CLUSTER_ALGORITHM = ['leiden', 'louvain']
sc.settings.set_figure_params(dpi=120)
sc.settings.verbosity = 0

logger = getlogger(__name__)
logger.setLevel(10)


def cluster(ctx, matrix, outdir, sample, barcodes, genes):
    sample_outdir = Path(outdir, sample, '06.cluster')
    sample_outdir.mkdir(parents=True, exist_ok=True)
    os.chdir(sample_outdir)

    adata_mtx = sc.read_mtx(matrix).T
    obs = pd.read_csv(barcodes, index_col=0, header=None)
    var = pd.read_csv(genes, index_col=0, header=None)
    obs.index.set_names(None, inplace=True)
    var.index.set_names(None, inplace=True)
    adata_mtx.obs = obs
    adata_mtx.var = var

    result_file = sample_outdir / f'{sample}.h5ad'

    adata_mtx.var_names_make_unique()
    [sc.pl.highest_expr_genes(adata_mtx, n_top=30, save=f'_{sample}.{filetype}') for filetype in FILETYPE]

    sc.pp.filter_cells(adata_mtx, min_genes=200)
    sc.pp.filter_genes(adata_mtx, min_cells=3)

    adata_mtx.var['mt'] = adata_mtx.var_names.str.upper().str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata_mtx, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    [sc.pl.violin(adata_mtx, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save=f'_{sample}.{filetype}') for filetype in FILETYPE]
    [sc.pl.scatter(adata_mtx, x='total_counts', y='pct_counts_mt', save=f'_pct_mt_{sample}.{filetype}') for filetype in FILETYPE]
    [sc.pl.scatter(adata_mtx, x='total_counts', y='n_genes_by_counts', save=f'_n_genes_{sample}.{filetype}') for filetype in FILETYPE]

    adata_mtx = adata_mtx[adata_mtx.obs.n_genes_by_counts < 2500, :]
    adata_mtx = adata_mtx[adata_mtx.obs.pct_counts_mt < 5, :]
    sc.pp.normalize_total(adata_mtx, target_sum=1e4)
    sc.pp.log1p(adata_mtx)
    sc.pp.highly_variable_genes(adata_mtx, min_mean=0.0125, max_mean=3, min_disp=0.5)

    [sc.pl.highly_variable_genes(adata_mtx, save=f'_{sample}.{filetype}') for filetype in FILETYPE]

    adata_mtx.raw = adata_mtx
    adata_mtx = adata_mtx[:, adata_mtx.var.highly_variable]

    sc.pp.regress_out(adata_mtx, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata_mtx, max_value=10)
    sc.tl.pca(adata_mtx, svd_solver='arpack')

    sc.pp.neighbors(adata_mtx, n_neighbors=20, n_pcs=40)
    sc.tl.umap(adata_mtx)
    sc.tl.leiden(adata_mtx)
    sc.tl.louvain(adata_mtx)
    [sc.pl.umap(adata_mtx, color=algo, save=f'_{algo}_{sample}.{filetype}') for filetype in FILETYPE for algo in CLUSTER_ALGORITHM]

    adata_mtx.write(result_file, compression='gzip')

    stat_info = []
    img = []
    with open(sample_outdir / 'stat.json', mode='w', encoding='utf-8') as f:
        pngs = Path(sample_outdir / 'figures').rglob('*.png')
        for png in pngs:
            img.append(
                {
                    'path': str(png.resolve()),
                    'name': png.name
                }
            )
        json.dump(stat_info, f)

    # report
    logger.info('generate report start!')
    Reporter(name='cluster', stat_file=sample_outdir / 'stat.json', outdir=sample_outdir.parent, img=img)

    logger.info('generate report done!')
