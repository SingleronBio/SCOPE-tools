# -*- coding: utf-8 -*-

import csv
from collections import defaultdict
from itertools import groupby
from ._count import umi_reads_downsample as downsample
import matplotlib.pyplot as plt
import pandas as pd
import pysam
from scipy.io import mmwrite
from scipy.sparse import coo_matrix

from .report import Reporter
from .utils import getlogger, cached_property

logger = getlogger(__name__)
logger.setLevel(10)


class CellGeneUmiSummary(object):

    def __init__(self, file, outdir, sample, cell_num=3000):
        """

        :param file: Path
        :param outdir: Path
        :param sample: str
        :param cell_num: int
        """

        self.file = file
        self.sample = sample
        self.cell_num = cell_num

        self.seq_df = pd.read_csv(self.file, index_col=[0])
        self.cell_df = None
        self.nth = None
        self.threshold = None
        self.valid_cell = None
        self.saturations = downsample(self.seq_df)
        self.count_info = {}
        self.umi_info = {}

        self.pdf = outdir / 'barcode_filter_magnitude.pdf'
        self.marked_counts_file = outdir / f'{self.sample}_counts.csv'
        self.matrix_file = outdir / f'{self.sample}_matrix.mtx'
        self.matrix_cellbarcode_file = outdir / f'{self.sample}_barcodes.tsv'
        self.matrix_gene_file = outdir / f'{self.sample}_genes.tsv'

        self.call_cells()
        self.plot_umi_cell()
        self.generate_matrix()
        self.generate_count_summary()
        self.generate_umi_summary()

    @cached_property
    def gene_cell_matrix(self):
        matrix = self.seq_df.loc[self.seq_df['mark'] > 0, :].pivot_table(index='geneID', columns='Barcode', aggfunc={'UMI': 'count'}).fillna(0).astype(int)
        matrix.columns = matrix.columns.droplevel(level=0)
        return matrix

    @cached_property
    def cell_describe(self):
        return self.cell_df.loc[self.cell_df['mark'] > 0, :].describe()

    @cached_property
    def cell_total_genes(self):
        return self.seq_df.loc[self.seq_df['mark'] > 0, 'geneID'].nunique()

    @cached_property
    def cell_reads_count(self):
        return self.seq_df.loc[self.seq_df['mark'] > 0, 'count'].sum()

    @cached_property
    def reads_mapped_to_transcriptome(self):
        return self.seq_df['count'].sum()

    def call_cells(self):
        self.cell_df = self.seq_df.pivot_table(
            index=['Barcode'],
            aggfunc={
                'geneID': 'nunique',
                'UMI': 'count',
                'count': ['sum', lambda x: sum(x[x > 1])]
            }
        )
        self.cell_df.columns = ['UMI', 'UMI2', 'read_count', 'geneID']

        self.nth = max(0, int(self.cell_num * 0.01) - 1)
        self.threshold = max(1, int(self.cell_df['UMI'].nlargest(self.nth)[-1] * 0.1))
        self.valid_cell = self.cell_df[self.cell_df['UMI'] > self.threshold].index

        self.seq_df.loc[:, 'mark'] = 0
        self.seq_df.loc[self.seq_df.index.isin(self.valid_cell), 'mark'] = 1

        self.cell_df.loc[:, 'mark'] = 0
        self.cell_df.loc[self.cell_df.index.isin(self.valid_cell), 'mark'] = 1

        self.cell_df.to_csv(self.marked_counts_file)

    def plot_umi_cell(self):
        plt.plot(self.cell_df['UMI'].sort_values(ascending=False))
        plt.hlines(self.threshold, 0, self.valid_cell.shape[0], linestyle='dashed')
        plt.vlines(self.valid_cell.shape[0], 0, self.threshold, linestyle='dashed')
        plt.xlabel('cell count')
        plt.ylabel('UMI num')
        plt.title(
            f'expected cell num:{self.cell_num}\nUMI threshold:{self.threshold}\ncell num:{self.valid_cell.shape[0]}')
        plt.loglog()
        plt.savefig(self.pdf)

    def generate_matrix(self):
        self.gene_cell_matrix.columns.to_series().to_csv(self.matrix_cellbarcode_file, index=False, header=False)
        self.gene_cell_matrix.index.to_series().to_csv(self.matrix_gene_file, index=False, header=False)
        mmwrite(str(self.matrix_file), coo_matrix(self.gene_cell_matrix))

    def generate_count_summary(self):
        attrs = [
            'Cells_number',
            'Saturation',
            'Mean_Reads',
            'Median_UMIs',
            'Total_Genes',
            'Median_Genes',
            'fraction_reads_in_cells',
            'mean_reads_per_cell',
        ]
        vals = [
            f"{int(self.cell_describe.loc['count', 'read_count']):d}",
            f"{self.saturations.loc[1, 'saturation']:.2%}",
            f"{int(self.cell_describe.loc['mean', 'read_count']):d}",
            f"{int(self.cell_describe.loc['50%', 'UMI']):d}",
            f"{int(self.cell_total_genes):d}",
            f"{int(self.cell_describe.loc['50%', 'geneID'])}",
            f"{self.cell_reads_count / self.reads_mapped_to_transcriptome:.2f}",
            f"{int(self.reads_mapped_to_transcriptome / self.cell_describe.loc['count', 'read_count']):d}"
        ]
        for attr, val in zip(attrs, vals):
            self.count_info[attr] = val

    def generate_umi_summary(self):
        attrs = [
            'percentile',
            'MedianGeneNum',
            'Saturation',
            'CB_num',
            'Cells',
            'UB_num',
            'Background'
        ]
        vals = [
            self.saturations.index.astype(float).to_list(),
            list(map(int, self.saturations.iloc[:, 0].astype(float).to_list())),
            (100 * self.saturations.iloc[:, 1].astype(float)).to_list(),
            self.cell_df[self.cell_df['mark'] > 0].shape[0],
            self.cell_df[self.cell_df['mark'] > 0]['UMI'].sort_values(ascending=False).to_list(),
            self.cell_df[self.cell_df['mark'] < 1].shape[0],
            self.cell_df[self.cell_df['mark'] < 1]['UMI'].sort_values(ascending=False).to_list(),
        ]

        for attr, val in zip(attrs, vals):
            self.umi_info[attr] = val


class Cell(object):

    @staticmethod
    def is_distance(x, y):
        return sum([True for i, j in zip(x, y) if i != j]) == 1

    def __init__(self, cell, seqs):
        self.cell = cell
        self.seqs = seqs
        self.gene_umi = defaultdict(lambda: defaultdict(int))
        self.parse_seqs()
        self.correct_umi()

    def parse_seqs(self):
        for seq in self.seqs:
            if seq.has_tag('XT'):
                gene_id = seq.get_tag('XT')
                cell, umi = seq.qname.split('_')[:2]
                self.gene_umi[gene_id][umi] += 1

    def correct_umi(self, percent=0.1):
        for gene in self.gene_umi:
            umis = self.gene_umi[gene]
            _umis = umis.copy()
            # only one umi on the gene,do not correct
            if len(umis) == 1:
                continue

            umis_sorted = sorted(umis.keys(), key=lambda x: umis[x], reverse=True)
            mis_umi = []

            for i in range(len(umis) - 1, 0, -1):
                umi_min = umis_sorted[i]
                if umi_min in mis_umi:
                    continue
                for j in range(i):
                    umi = umis_sorted[j]
                    if umis[umi_min] / umis[umi] > percent:
                        continue
                    if self.is_distance(umi_min, umi):
                        mis_umi.append(umi_min)
                        _umis[umi] += umis[umi_min]

            self.gene_umi[gene] = _umis


def count(ctx, bam, sample, outdir, cells):
    sample_outdir = outdir / sample / '05.count'
    sample_outdir.mkdir(parents=True, exist_ok=True)

    # umi correction
    logger.info('UMI count start!')
    count_detail_file = sample_outdir / f'{sample}_count_detail.csv'
    with pysam.AlignmentFile(bam, mode='r') as p, open(count_detail_file, mode='w') as f:
        headers = ['Barcode', 'geneID', 'UMI', 'count']

        f_csv = csv.writer(f)
        f_csv.writerow(headers)

        for cell, seqs in groupby(p, lambda seq: seq.qname.split('_', 1)[0]):
            c = Cell(cell, seqs)
            if c.gene_umi:
                for gene in c.gene_umi:
                    for umi in c.gene_umi[gene]:
                        f_csv.writerow((c.cell, gene, umi, c.gene_umi[gene][umi]))
    logger.info('UMI count done!')

    # call cells
    logger.info('cell count start!')
    cell_gene = CellGeneUmiSummary(file=count_detail_file, outdir=sample_outdir, sample=sample, cell_num=cells)
    logger.info('cell count done!')

    # report
    logger.info('generate report start!')
    Reporter(name='count', stat_json=cell_gene.count_info, outdir=sample_outdir.parent)
    Reporter(name='UMI', stat_json=cell_gene.umi_info, outdir=sample_outdir.parent)
    logger.info('generate report done!')
