# -*- coding: utf-8 -*-

from pathlib import Path
from scopetools.utils import getlogger
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

logger = getlogger(__name__)
logger.setLevel(10)


class CellGeneUmiSummary(object):

    def __init__(self, file: Path, outdir: Path, sample: str, cell_num: int = 3000):
        self.seq_df = pd.read_csv(file)
        self.cell_num = cell_num
        self.sample = sample

        self.pdf = outdir / 'barcode_filter_magnitude.pdf'
        self.marked_counts_file = outdir / f'{self.sample}_counts.csv'
        self.matrix_file = outdir / f'{self.sample}_matrix.xls'
        self.matrix_cellbarcode_file = outdir / f'{self.sample}_cellbarcode.tsv'
        self.matrix_gene_file = outdir / f'{self.sample}_genes.tsv'

        self.threshold, self.valid_cell, self.cell_df = self.call_cells()
        self.plot_umi_cell()
        self.cellbarcode_describe, self.cellbarcode_total_genes, self.cellbarcode_reads_count, self.reads_mapped_to_transcriptome = self.gen_matrix()
        self.save()

    def call_cells(self, col='UMI'):
        cell_df = self.seq_df.groupby('Barcode').agg(
            {
                'count': ['sum', lambda x: sum(x[x > 1])],
                'UMI': 'count',
                'geneID': 'nunique'
            }
        )
        cell_df.columns = ['read_count', 'UMI2', 'UMI', 'geneID']
        cell_df.sort_values(by=col, ascending=False, inplace=True)
        nth = max(0, int(self.cell_num * 0.01) - 1)
        threshold = max(1, int(cell_df.iloc[nth][col] * 0.1))
        valid_cell = cell_df[cell_df[col] > threshold].index
        cell_df.loc[:, 'mark'] = 'UB'
        cell_df.loc[cell_df.index.isin(valid_cell), 'mark'] = 'CB'
        cell_df.to_csv(self.marked_counts_file)
        return threshold, valid_cell, cell_df

    def plot_umi_cell(self, col='UMI'):
        plt.plot(self.cell_df[col])
        plt.hlines(self.threshold, 0, self.valid_cell.shape[0], linestyle='dashed')
        plt.vlines(self.valid_cell.shape[0], 0, self.threshold, linestyle='dashed')
        plt.xlabel('cell count')
        plt.ylabel('UMI num')
        plt.title(f'expected cell num:{self.cell_num}\nUMI threshold:{self.threshold}\ncell num:{self.valid_cell.shape[0]}')
        plt.loglog()
        plt.savefig(self.pdf)

    def gen_matrix(self):
        self.seq_df.loc[:, 'mark'] = 'UB'
        self.seq_df.loc[self.seq_df['Barcode'].isin(self.valid_cell), 'mark'] = 'CB'

        cellbarcode_describe = self.cell_df.loc[self.cell_df['mark'] == 'CB', :].describe()
        cellbarcode_total_genes = self.seq_df.loc[self.seq_df['mark'] == 'CB', 'geneID'].nunique()
        cellbarcode_reads_count = self.seq_df.loc[self.seq_df['mark'] == 'CB', 'count'].sum()
        reads_mapped_to_transcriptome = self.seq_df['count'].sum()

        table = self.cell_df.loc[self.cell_df['mark'] == 'CB', :].pivot_table(index='geneID', columns='Barcode', values='UMI', aggfunc=len).fillna(0).astype(int)
        table.fillna(0).to_csv(self.matrix_file, sep='\t')
        table.columns.to_series().to_csv(self.matrix_cellbarcode_file, index=False, sep='\t')
        table.index.to_series().to_csv(self.matrix_gene_file, index=False, sep='\t')
        mmwrite(str(self.matrix_file), csr_matrix(table.fillna(0)))
        return cellbarcode_describe, cellbarcode_total_genes, cellbarcode_reads_count, reads_mapped_to_transcriptome

    def save(self):
        self.seq_df.to_csv(self.marked_counts_file)


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
    sample_outdir = Path(outdir, sample, '05.count')
    sample_outdir.mkdir(parents=True, exist_ok=True)

    # umi correction
    logger.info('UMI count start!')
    count_detail_file = sample_outdir / f'{sample}_count_detail.csv'
    # with pysam.AlignmentFile(bam, mode='r') as p, open(count_detail_file, mode='w') as f:
    #     headers = ['Barcode', 'geneID', 'UMI', 'count']
    #
    #     f_csv = csv.writer(f)
    #     f_csv.writerow(headers)
    #
    #     for cell, seqs in groupby(p, lambda seq: seq.qname.split('_', 1)[0]):
    #         c = Cell(cell, seqs)
    #         if c.gene_umi:
    #             for gene in c.gene_umi:
    #                 for umi in c.gene_umi[gene]:
    #                     f_csv.writerow((c.cell, gene, umi, c.gene_umi[gene][umi]))
    # logger.info('UMI count done!')

    # call cells
    pdf = sample_outdir / 'barcode_filter_magnitude.pdf'
    marked_counts_file = sample_outdir / f'{sample}_counts.txt'
    umi_plot = CellGeneUmiSummary(file=count_detail_file, outdir=sample_outdir, sample=sample, cell_num=cells)
