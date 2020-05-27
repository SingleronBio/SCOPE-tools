# -*- coding: utf-8 -*-
from pathlib import Path

import click

from scopetools.utils import BarcodeType, AdapterType, MultipleOption, MutuallyExclusiveOption

__version__ = '0.1.0'
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(chain=True, context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__, prog_name='scope-tools')
def cli():
    """
    Single Cell Omics Preparation Entity Tools
    """
    pass


@cli.command(name='barcode', short_help="extract barcode and umi short help")
@click.option('--fq1', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True, help="read1 fq file")
@click.option('--fq2', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True, help="read2 fq file")
@click.option('--sample', type=click.STRING, required=True, help="sample name")
@click.option('--outdir', type=click.Path(file_okay=False, dir_okay=True, writable=True), required=True, help="output dir")
@click.option('--bctype', cls=MutuallyExclusiveOption, type=click.Choice(['SCOPEv2', 'SCOPEv1', '10X', 'Dropseq', 'inDrop', 'BD', 'other']), mutually_exclusive=['pattern', 'whitelist', 'linkers'], default=None, help='bctype help')
@click.option('--pattern', cls=MutuallyExclusiveOption, type=BarcodeType(), mutually_exclusive=['bctype'], default=None, help="read1 pattern, C: cellbarcode, L: linker, U: UMI, T: polyT")
@click.option('--whitelist', cls=MutuallyExclusiveOption, type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), mutually_exclusive=['bctype'], default=None, help="cell barcode list")
@click.option('--linkers', multiple=True, cls=MutuallyExclusiveOption, type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), mutually_exclusive=['bctype'], default=None, help="linkers")
@click.option('--lowQual', type=click.IntRange(0, 30), default=14, show_default=True, help="max phred of base as low quality")
@click.option('--lowNum', type=click.INT, default=2, show_default=True, help="max number with low quality allowed")
@click.pass_context
def barcode_pipe(ctx, fq1, fq2, sample, outdir, bctype, pattern, whitelist, linkers, lowqual, lownum):
    """
    extract barcode and umi description
    """
    if bctype:
        if bctype == 'SCOPEv2':
            pattern = 'C8L16C8L16C8U8T18'
            whitelist = Path(__file__).resolve().parent.joinpath(f'extra/whitelist/SCOPEv2/bclist')
            linkers = [Path(__file__).resolve().parent.joinpath(f'extra/whitelist/SCOPEv2/{linker}') for linker in ['linker1', 'linker2']]
        elif bctype == 'SCOPEv1':
            pattern = 'C12U8T30'
            whitelist = None
            linkers = None
        elif bctype == '10X':
            pattern = 'C16U12T30'
            whitelist = None
            linkers = None
    elif pattern:
        l_num = pattern.count('L')
        if l_num:
            if len(linkers) != pattern.count('L'):
                raise click.BadParameter("{} is not a valid adapter pattern".format(pattern))
    else:
        raise click.BadParameter("Error: Illegal usage: [bctype] or [pattern whitelist linkers] must have one")
    from scopetools.barcode import barcode
    barcode(ctx=ctx, fq1=fq1, fq2=fq2, sample=sample, outdir=outdir, pattern=pattern, whitelist=whitelist, linkers=linkers, lowqual=lowqual, lownum=lownum)


@cli.command(name='cutadapt', short_help="cutadapt short help")
@click.option('--fq', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True, help="fq help")
@click.option('--sample', type=click.STRING, required=True, help="sample help")
@click.option('--outdir', type=click.Path(file_okay=False, dir_okay=True, writable=True), required=True, help="outdir help")
@click.option('--adapter', cls=MultipleOption, type=AdapterType(), nargs=-1, default=['polyT=A{18}', 'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'], help="adapter help")
@click.option('--minimum-length', type=click.INT, default=20, show_default=True, help="minimum length help")
@click.option('--nextseq-trim', type=click.INT, default=20, show_default=True, help="nextseq trim help")
@click.option('--overlap', type=click.INT, default=5, show_default=True, help="overlap help")
@click.option('--thread', type=click.INT, default=2, show_default=True, help="thread help")
@click.pass_context
def cutadapt_pipe(ctx, fq, sample, outdir, adapter, minimum_length, nextseq_trim, overlap, thread):
    """
    cutadapt description
    """
    from scopetools.cutadapt import cutadapt
    cutadapt(ctx, fq, sample, outdir, adapter, minimum_length, nextseq_trim, overlap, thread)


@cli.command(name="STAR", short_help="STAR short help")
@click.option('--fq', type=click.Path(exists=True, file_okay=True, readable=True), required=True, help="fq help")
@click.option('--refFlat', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True, default='/SGRNJ/Public/Database/genome/homo_mus/homo_mus.refFlat', help="refFlat help")
@click.option('--genomeDir', type=click.Path(file_okay=False, dir_okay=True, readable=True), required=True, default='/SGRNJ/Public/Database/genome/homo_mus', help="genome help")
@click.option('--sample', type=click.STRING, required=True, help="sample help")
@click.option('--outdir', type=click.Path(file_okay=False, dir_okay=True, writable=True), required=True, help="outdir help")
@click.option('--readFilesCommand', type=click.Choice(['zcat', 'bzcat', '-']), default='zcat', show_default=True, help='readFilesCommand help')
@click.option('--runThreadN', type=click.INT, default=2, show_default=True, help='runThreadN help')
@click.pass_context
def star_pipe(ctx, fq, refflat, genomedir, sample, outdir, readfilescommand, runthreadn):
    """
    STAR description
    """
    from scopetools.star import star
    star(ctx, fq, refflat, genomedir, sample, outdir, readfilescommand, runthreadn)


@cli.command(name='featureCounts', short_help="featureCounts short help")
@click.option('--input', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True, help="input help")
@click.option('--annot', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True, help="annot help")
@click.option('--sample', type=click.STRING, required=True, help="sample help")
@click.option('--outdir', type=click.Path(file_okay=False, dir_okay=True, writable=True), required=True, help="outdir help")
@click.option('--format', type=click.Choice(['SAM', 'BAM', 'CORE']), default='BAM', show_default=True, help="format help")
@click.option('--nthreads', type=click.IntRange(1, 32), default=2, show_default=True, help="nthreads help")
@click.pass_context
def featurecounts_pipe(ctx, input, annot, sample, outdir, format, nthreads):
    """
    featureCounts description
    """
    from scopetools.featurecounts import featurecounts
    featurecounts(ctx, input, annot, sample, outdir, format, nthreads)


@cli.command(name='count', short_help="count short help")
@click.option('--bam', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True, help="bam help")
@click.option('--sample', type=click.STRING, required=True, help="sample help")
@click.option('--outdir', type=click.Path(file_okay=False, dir_okay=True, writable=True), required=True, help="outdir help")
@click.option('--cells', type=click.INT, default=3000, show_default=True, help="sample help")
@click.pass_context
def count_pipe(ctx, bam, sample, outdir, cells):
    """
    count description
    """
    from scopetools.count import count
    count(ctx, bam, sample, outdir, cells)


@cli.command(name='cluster', short_help="cluster short help")
@click.option('--outdir', type=click.Path(file_okay=False, dir_okay=True, writable=True), required=True, help="outdir help")
@click.option('--sample', type=click.STRING, required=True, help="sample help")
@click.option('--matrix', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True, help="matrix help")
@click.option('--barcodes', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True, help="barcodes help")
@click.option('--genes', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True, help="genes help")
@click.option('--n_top', type=click.INT, default=30, show_default=True, help="Number of top")
@click.option('--min_genes', type=click.INT, default=200, show_default=True, help="Minimum number of genes expressed required for a cell to pass filtering")
@click.option('--min_cells', type=click.INT, default=3, show_default=True, help="Minimum number of cells expressed required for a gene to pass filtering")
@click.option('--n_genes_by_counts', type=click.INT, default=2500, show_default=True, help="Minimum number of expressed genes required for a cell to pass filtering")
@click.option('--pct_counts_mt', type=click.INT, default=5, show_default=True, help="Maximum pct_counts_mt required for a cell to pass filtering")
@click.option('--exclude_highly_expressed', type=click.BOOL, default=False, show_default=True, help="Exclude (very) highly expressed genes for the computation of the normalization factor (size factor) for each cell")
@click.option('--max_fraction', type=click.FLOAT, default=0.05, show_default=True, help="Consider cells as highly expressed that have more counts than max_fraction of the original total counts in at least one cell")
@click.option('--n_top_genes', type=click.INT, default=None, show_default=True, help="Number of highly-variable genes to keep.")
@click.option('--max_value', type=click.FLOAT, default=None, show_default=True, help="Clip (truncate) to this value after scaling")
@click.option('--n_neighbors', type=click.IntRange(2, 100), default=15, show_default=True, help="The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation")
@click.option('--n_pcs', type=click.INT, default=None, show_default=True, help="Use this many PCs")
@click.pass_context
def cluster_pipe(ctx, matrix, outdir, sample, barcodes, genes, n_top, min_genes, min_cells, n_genes_by_counts, pct_counts_mt, exclude_highly_expressed, max_fraction, n_top_genes, max_value, n_neighbors, n_pcs):
    """
    cluster description
    """
    click.echo('cluster pipeline')
    from scopetools.cluster import cluster
    cluster(ctx, matrix, outdir, sample, barcodes, genes, n_top, min_genes, min_cells, n_genes_by_counts, pct_counts_mt, exclude_highly_expressed, max_fraction, n_top_genes, max_value, n_neighbors, n_pcs)


@cli.command(name='run', help="run short help")
@click.option('--fq1', type=click.Path(exists=True, file_okay=True, readable=True), required=True, help="[barcode], read1 fq file")
@click.option('--fq2', type=click.Path(exists=True, file_okay=True, readable=True), required=True, help="[barcode], read2 fq file")
@click.option('--sample', type=click.STRING, required=True, help="[barcode], sample name")
@click.option('--outdir', type=click.Path(dir_okay=True, writable=True), required=True, help="[barcode], output dir")
@click.option('--pattern', type=BarcodeType(), default='C6L15C6L15C6U6T30', show_default=True, help="[barcode], read1 pattern, C: cellbarcode, L: linker, U: UMI, T: polyT")
@click.option('--whitelist', type=click.Path(exists=True, file_okay=True, readable=True), default=Path(__file__).resolve().parent.joinpath('extra/whitelist/scope/bclist'), show_default=False, help="[barcode], cell barcode list")
@click.option('--linkers', multiple=True, type=click.Path(exists=True, file_okay=True, readable=True), default=[Path(__file__).resolve().parent.joinpath('extra/whitelist/scope/linker1'), Path(__file__).resolve().parent.joinpath('extra/whitelist/scope/linker2')], show_default=False, help="[barcode], linkers")
@click.option('--lowQual', type=click.IntRange(0, 30), default=14, show_default=True, help="[barcode], max phred of base as low quality")
@click.option('--lowNum', type=click.INT, default=2, show_default=True, help="[barcode], max number with low quality allowed")
@click.option('--adapter', cls=MultipleOption, type=AdapterType(), nargs=-1, default=['polyT=A{18}', 'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'], help="[cutadapt], adapter help")
@click.option('--minimum-length', type=click.INT, default=20, show_default=True, help="[cutadapt], minimum length help")
@click.option('--nextseq-trim', type=click.INT, default=20, show_default=True, help="[cutadapt], nextseq trim help")
@click.option('--overlap', type=click.INT, default=5, show_default=True, help="[cutadapt], overlap help")
@click.option('--refFlat', type=click.Path(exists=True, file_okay=True, readable=True), required=True, help="[STAR], refFlat help")
@click.option('--genomeDir', type=click.Path(dir_okay=True, readable=True), required=True, help="[STAR], genome help")
@click.option('--readFilesCommand', type=click.Choice(['zcat', 'bzcat', '-']), default='zcat', show_default=True, help='[STAR], readFilesCommand help')
@click.option('--runThreadN', type=click.INT, default=2, show_default=True, help='[STAR], runThreadN help')
@click.option('--annot', type=click.Path(exists=True, file_okay=True, readable=True), required=True, help="[featureCounts], annot help")
@click.option('--format', type=click.Choice(['SAM', 'BAM', 'CORE']), default='BAM', show_default=True, help="[featureCounts], format help")
@click.option('--nthreads', type=click.IntRange(1, 32), default=2, show_default=True, help="[featureCounts], nthreads help")
@click.option('--cells', type=click.INT, default=3000, show_default=True, help="[count], sample help")
@click.pass_context
def run_pipe(ctx, fq1, fq2, sample, outdir, pattern, whitelist, linkers, lowqual, lownum, adapter, minimum_length, nextseq_trim, overlap, refflat, genomedir, readfilescommand, runthreadn, annot, format, nthreads, cells):
    """
    run description
    """
    click.echo('run pipeline')
    ctx.invoke(barcode_pipe)
    ctx.invoke(cutadapt_pipe)
    ctx.invoke(star_pipe)
    ctx.invoke(featurecounts_pipe)
    ctx.invoke(count_pipe)


if __name__ == '__main__':
    cli()
