# -*- coding: utf-8 -*-
import pandas as pd


def umi_reads_downsample(seq_df):
    """

    :param seq_df:
    :return: saturations
    """
    saturations = pd.DataFrame(columns=['percent', 'median_gene_num', 'saturation']).set_index('percent')
    all_seq_df = seq_df.reset_index().set_index(['Barcode', 'geneID', 'UMI', 'mark']).index.repeat(seq_df['count']).to_frame().set_index(['Barcode'])
    saturations.loc[0, :] = [0, 0]
    for i in range(1, 11):
        sample_df = all_seq_df.sample(frac=i / 10)
        sample_df = sample_df.loc[sample_df['mark'] > 0]
        total = sample_df['UMI'].count()
        gene_num_median = sample_df.pivot_table(index='Barcode', aggfunc={'geneID': 'nunique'})['geneID'].median()
        sample_df = sample_df.pivot_table(index=['Barcode', 'geneID', 'UMI'], aggfunc={'UMI': 'count'})
        repeat = sample_df.loc[sample_df['UMI'] > 1, 'UMI'].sum()
        saturation = repeat / total
        saturations.loc[i / 10, :] = [gene_num_median, saturation]
    return saturations


def umi_count_downsample(seq_df):
    """

    :param seq_df:
    :return: saturations
    """
    saturations = pd.DataFrame(columns=['percent', 'median_gene_num', 'saturation']).set_index('percent')
    all_seq_df = seq_df.reset_index().set_index(['Barcode', 'geneID', 'UMI', 'mark']).index.repeat(seq_df['count']).to_frame().set_index(['Barcode'])
    saturations.loc[0, :] = [0, 0]
    for i in range(1, 11):
        sample_df = all_seq_df.sample(frac=i / 10)
        sample_df = sample_df.loc[sample_df['mark'] > 0]
        tmp = sample_df.pivot_table(index=['Barcode', 'geneID', 'UMI'], aggfunc={'mark': 'count'}).reset_index().set_index(['Barcode'])
        total = tmp.pivot_table(index=['Barcode'], aggfunc={'mark': 'count'})['mark'].sum()
        repeat = tmp[tmp['mark'] > 1].pivot_table(index=['Barcode'], aggfunc={'mark': 'count'})['mark'].sum()
        sample_df_pivot = sample_df.pivot_table(index=['Barcode'], aggfunc={'UMI': 'count', 'geneID': 'nunique'})
        gene_num_median = sample_df_pivot['geneID'].median()
        saturation = repeat / total
        saturations.loc[i / 10, :] = [gene_num_median, saturation]
    return saturations
