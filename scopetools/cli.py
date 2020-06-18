# -*- coding: utf-8 -*-
from pathlib import Path

import click

from ._version import get_versions
from .utils import BarcodeType, AdapterType, MutuallyExclusiveOption, str2path

__version__ = get_versions()['version']
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(chain=True, context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__, prog_name='SCOPE-tools')
def cli():
    """
    Single Cell Omics Preparation Entity Tools
    """
    pass


def sample_param(func):
    func = click.option('--sample', type=click.STRING, required=True, help="sample name")(func)
    func = click.option('--description', type=click.STRING, default='scRNA-Seq scope', show_default=True, help='description help')(func)
    func = click.option('--version', type=click.STRING, default=__version__, show_default=True, help='version help')(func)
    func = click.option('--outdir', type=click.Path(file_okay=False, dir_okay=True, writable=True), callback=str2path, required=True, help="output dir")(func)
    return func


def barcode_param(func):
    func = click.option('--fq1', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), multiple=True, callback=str2path, required=True, help="read1 fq file")(func)
    func = click.option('--fq2', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), multiple=True, callback=str2path, required=True, help="read2 fq file")(func)
    func = click.option('--bctype', cls=MutuallyExclusiveOption, type=click.Choice(['SCOPEv2', 'SCOPEv1', '10X', 'Dropseq', 'inDrop', 'BD', 'other']), mutually_exclusive=['pattern', 'whitelist', 'linkers'], default=None, help='bctype help')(func)
    func = click.option('--pattern', cls=MutuallyExclusiveOption, type=BarcodeType(), mutually_exclusive=['bctype'], default=None, help="read1 pattern, C: cellbarcode, L: linker, U: UMI, T: polyT")(func)
    func = click.option('--whitelist', cls=MutuallyExclusiveOption, type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), callback=str2path, mutually_exclusive=['bctype'], default=None, help="cell barcode list")(func)
    func = click.option('--linker', cls=MutuallyExclusiveOption, type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), callback=str2path, mutually_exclusive=['bctype'], default=None, help="linkers")(func)
    func = click.option('--lowQual', type=click.IntRange(0, 30), default=14, show_default=True, help="max phred of base as low quality")(func)
    func = click.option('--lowNum', type=click.INT, default=2, show_default=True, help="max number with low quality allowed")(func)
    return func


def cutadapt_param(func):
    func = click.option('--adapter', type=AdapterType(), default=['A{18}'], multiple=True, help="adapter help")(func)
    func = click.option('--minimum-length', type=click.INT, default=20, show_default=True, help="minimum length help")(func)
    func = click.option('--nextseq-trim', type=click.INT, default=20, show_default=True, help="nextseq trim help")(func)
    func = click.option('--overlap', type=click.INT, default=5, show_default=True, help="overlap help")(func)
    func = click.option('--thread', type=click.INT, default=2, show_default=True, help="thread help")(func)
    return func


def star_param(func):
    func = click.option('--refFlat', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), callback=str2path, required=True, help="refFlat help")(func)
    func = click.option('--genomeDir', type=click.Path(file_okay=False, dir_okay=True, readable=True), callback=str2path, required=True, help="genome help")(func)
    func = click.option('--readFilesCommand', type=click.Choice(['zcat', 'bzcat', '-']), default='zcat', show_default=True, help='readFilesCommand help')(func)
    func = click.option('--runThreadN', type=click.INT, default=2, show_default=True, help='runThreadN help')(func)
    return func


def featurecounts_param(func):
    func = click.option('--annot', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), callback=str2path, required=True, help="annot help")(func)
    func = click.option('--format', type=click.Choice(['SAM', 'BAM', 'CORE']), default='BAM', show_default=True, help="format help")(func)
    func = click.option('--nthreads', type=click.IntRange(1, 32), default=2, show_default=True, help="nthreads help")(func)
    func = click.option('--type', type=click.STRING, default='exon', show_default=True, help="specify feature type in GTF annotation")(func)
    return func


def count_param(func):
    func = click.option('--cells', type=click.INT, default=3000, show_default=True, help="sample help")(func)
    return func


def cluster_param(func):
    func = click.option('--n_top', type=click.INT, default=30, show_default=True, help="Number of top")(func)
    func = click.option('--min_genes', type=click.INT, default=200, show_default=True, help="Minimum number of genes expressed required for a cell to pass filtering")(func)
    func = click.option('--min_cells', type=click.INT, default=3, show_default=True, help="Minimum number of cells expressed required for a gene to pass filtering")(func)
    func = click.option('--n_genes_by_counts', type=click.INT, default=2500, show_default=True, help="Minimum number of expressed genes required for a cell to pass filtering")(func)
    func = click.option('--pct_counts_mt', type=click.INT, default=5, show_default=True, help="Maximum pct_counts_mt required for a cell to pass filtering")(func)
    func = click.option('--exclude_highly_expressed', type=click.BOOL, default=False, show_default=True, help="Exclude (very) highly expressed genes for the computation of the normalization factor (size factor) for each cell")(func)
    func = click.option('--max_fraction', type=click.FLOAT, default=0.05, show_default=True, help="Consider cells as highly expressed that have more counts than max_fraction of the original total counts in at least one cell")(func)
    func = click.option('--n_top_genes', type=click.INT, default=None, show_default=True, help="Number of highly-variable genes to keep.")(func)
    func = click.option('--max_value', type=click.FLOAT, default=None, show_default=True, help="Clip (truncate) to this value after scaling")(func)
    func = click.option('--n_neighbors', type=click.IntRange(2, 100), default=15, show_default=True, help="The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation")(func)
    func = click.option('--n_pcs', type=click.INT, default=None, show_default=True, help="Use this many PCs")(func)
    return func


@cli.command(name='sample', short_help="sample short help")
@sample_param
@click.option('--transcriptome', type=click.STRING, default='Unspecified', required=True, help="sample help")
@click.pass_context
def sample_pipe(ctx, sample, transcriptome, description, version, outdir, *args, **kwargs):
    """
    sample description
    """
    from .sample_description import sample_description
    sample_description(ctx, sample, transcriptome, description, version, outdir)


@cli.command(name='barcode', short_help="extract barcode and umi short help")
@barcode_param
@click.option('--sample', type=click.STRING, required=True, help="sample help")
@click.option('--outdir', type=click.Path(file_okay=False, dir_okay=True, writable=True), callback=str2path, required=True, help="outdir help")
@click.pass_context
def barcode_pipe(ctx, fq1, fq2, sample, outdir, bctype, pattern, whitelist, linker, lowqual, lownum, *args, **kwargs):
    """
    extract barcode and umi description
    """
    if bctype:
        if bctype == 'SCOPEv2':
            pattern = 'C8L16C8L16C8L1U8T18'
            whitelist = Path(__file__).resolve().parent.joinpath(f'extra/whitelist/SCOPEv2/bclist')
            linker = Path(__file__).resolve().parent.joinpath(f'extra/whitelist/SCOPEv2/linker')
        elif bctype == 'SCOPEv1':
            pattern = 'C12U8T30'
            whitelist = None
            linker = None
        elif bctype == '10X':
            pattern = 'C16U12T30'
            whitelist = None
            linker = None
    elif pattern:
        pass
        # raise click.BadParameter("{} is not a valid adapter pattern".format(pattern))
    else:
        raise click.BadParameter("Error: Illegal usage: [bctype] or [pattern whitelist linkers] must have one")
    from .barcode import barcode
    barcode(ctx=ctx, fq1s=fq1, fq2s=fq2, sample=sample, outdir=outdir, bctype=bctype, pattern=pattern, whitelist=whitelist, linker=linker, lowqual=lowqual, lownum=lownum)


@cli.command(name='cutadapt', short_help="cutadapt short help")
@cutadapt_param
@click.option('--fq', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), callback=str2path, required=True, help="fq help")
@click.option('--sample', type=click.STRING, required=True, help="sample help")
@click.option('--outdir', type=click.Path(file_okay=False, dir_okay=True, writable=True), callback=str2path, required=True, help="outdir help")
@click.pass_context
def cutadapt_pipe(ctx, fq, sample, outdir, adapter, minimum_length, nextseq_trim, overlap, thread, *args, **kwargs):
    """
    cutadapt description
    """
    from .cutadapt import cutadapt
    cutadapt(ctx, fq, sample, outdir, adapter, minimum_length, nextseq_trim, overlap, thread)


@cli.command(name="STAR", short_help="STAR short help")
@star_param
@click.option('--fq', type=click.Path(exists=True, file_okay=True, readable=True), callback=str2path, required=True, help="fq help")
@click.option('--outdir', type=click.Path(file_okay=False, dir_okay=True, writable=True), callback=str2path, required=True, help="outdir help")
@click.option('--sample', type=click.STRING, required=True, help="sample help")
@click.pass_context
def star_pipe(ctx, fq, refflat, genomedir, sample, outdir, readfilescommand, runthreadn, *args, **kwargs):
    """
    STAR description
    """
    from .star import star
    star(ctx, fq, refflat, genomedir, sample, outdir, readfilescommand, runthreadn)


@cli.command(name='featureCounts', short_help="featureCounts short help")
@click.option('--input', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True, help="input help")
@click.option('--annot', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True, help="annot help")
@click.option('--sample', type=click.STRING, required=True, help="sample help")
@click.option('--outdir', type=click.Path(file_okay=False, dir_okay=True, writable=True), callback=str2path, required=True, help="outdir help")
@click.pass_context
def featurecounts_pipe(ctx, input, annot, sample, outdir, format, nthreads, type, *args, **kwargs):
    """
    featureCounts description
    """
    from .featurecounts import featurecounts
    featurecounts(ctx, input, annot, sample, outdir, format, nthreads, type)


@cli.command(name='count', short_help="count short help")
@count_param
@click.option('--bam', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), callback=str2path, required=True, help="bam help")
@click.option('--sample', type=click.STRING, required=True, help="sample help")
@click.option('--outdir', type=click.Path(file_okay=False, dir_okay=True, writable=True), required=True, help="outdir help")
@click.option('--cells', type=click.INT, default=3000, show_default=True, help="sample help")
@click.pass_context
def count_pipe(ctx, bam, sample, outdir, cells, *args, **kwargs):
    """
    count description
    """
    from .count import count
    count(ctx, bam, sample, outdir, cells)


@cli.command(name='cluster', short_help="cluster short help")
@cluster_param
@click.option('--outdir', type=click.Path(file_okay=False, dir_okay=True, writable=True), callback=str2path, required=True, help="outdir help")
@click.option('--sample', type=click.STRING, required=True, help="sample help")
@click.option('--matrix', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), callback=str2path, required=True, help="matrix help")
@click.option('--barcodes', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), callback=str2path, required=True, help="barcodes help")
@click.option('--genes', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), callback=str2path, required=True, help="genes help")
@click.pass_context
def cluster_pipe(ctx, matrix, outdir, sample, barcodes, genes, n_top, min_genes, min_cells, n_genes_by_counts, pct_counts_mt, exclude_highly_expressed, max_fraction, n_top_genes, max_value, n_neighbors, n_pcs, *args, **kwargs):
    """
    cluster description
    """
    click.echo('cluster pipeline')
    from .cluster import cluster
    cluster(ctx, matrix, outdir, sample, barcodes, genes, n_top, min_genes, min_cells, n_genes_by_counts, pct_counts_mt, exclude_highly_expressed, max_fraction, n_top_genes, max_value, n_neighbors, n_pcs)


@cli.command(name='run', help="run short help")
@sample_param
@barcode_param
@cutadapt_param
@star_param
@featurecounts_param
@count_param
@cluster_param
@click.pass_context
def run_pipe(ctx, description, version, fq1, fq2, sample, outdir, bctype, pattern, whitelist, linker, lowqual, lownum, adapter, minimum_length, nextseq_trim, overlap, thread, refflat, genomedir, readfilescommand, runthreadn, annot, format, nthreads, type, cells, n_top, min_genes, min_cells, n_genes_by_counts, pct_counts_mt, exclude_highly_expressed, max_fraction, n_top_genes, max_value, n_neighbors, n_pcs, *args, **kwargs):
    """
    run description
    """

    click.echo('run pipeline')

    # sample_pipe
    ctx.params['transcriptome'] = ctx.params['annot'].name.split('.')[0]
    ctx.invoke(sample_pipe, **ctx.params)

    # barcode_pipe
    ctx.invoke(barcode_pipe, **ctx.params)

    # cutadapt_pipe
    ctx.params['fq'] = ctx.params['outdir'] / ctx.params['sample'] / '01.barcode' / f'{ctx.params["sample"]}_2.fq.gz'
    if bctype == 'SCOPEv1':
        ctx.params['adapter'] = ['A{20}']
    ctx.invoke(cutadapt_pipe, **ctx.params)

    # star_pipe
    ctx.params['fq'] = ctx.params['outdir'] / ctx.params['sample'] / '02.cutadapt' / f'{ctx.params["sample"]}_clean_2.fq.gz'
    ctx.invoke(star_pipe, **ctx.params)

    # featurecounts_pipe
    ctx.params['input'] = ctx.params['outdir'] / ctx.params['sample'] / '03.STAR' / f'{ctx.params["sample"]}_Aligned.sortedByCoord.out.bam'
    ctx.invoke(featurecounts_pipe, **ctx.params)

    # count_pipe
    ctx.params['bam'] = ctx.params['outdir'] / ctx.params['sample'] / '04.featureCounts' / f'{ctx.params["sample"]}_name_sorted.bam'
    ctx.invoke(count_pipe, **ctx.params)

    # cluster_pipe
    ctx.params['matrix'] = ctx.params['outdir'] / ctx.params['sample'] / '05.count' / f'{ctx.params["sample"]}_matrix.mtx'
    ctx.params['barcodes'] = ctx.params['outdir'] / ctx.params['sample'] / '05.count' / f'{ctx.params["sample"]}_barcodes.tsv'
    ctx.params['genes'] = ctx.params['outdir'] / ctx.params['sample'] / '05.count' / f'{ctx.params["sample"]}_genes.tsv'
    ctx.invoke(cluster_pipe, **ctx.params)


if __name__ == '__main__':
    cli()
