# -*- coding: utf-8 -*-
import click
from scopetools.utils import BarcodeType, AdapterType, MultipleOption
from pathlib import Path

__version__ = '0.1.0'
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(chain=True, context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__, prog_name='scopetools')
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
@click.option('--pattern', type=BarcodeType(), default='C8L16C8L16C8U8T18', show_default=True, help="read1 pattern, C: cellbarcode, L: linker, U: UMI, T: polyT")
@click.option('--whitelist', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), default=Path(__file__).resolve().parent.joinpath('extra/whitelist/scope/bclist'), show_default=True, help="cell barcode list")
@click.option('--linkers', multiple=True, type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), default=[Path(__file__).resolve().parent.joinpath('extra/whitelist/scope/linker1'), Path(__file__).resolve().parent.joinpath('extra/whitelist/scope/linker2')], show_default=True, help="linkers")
@click.option('--lowQual', type=click.IntRange(0, 30), default=14, show_default=True, help="max phred of base as low quality")
@click.option('--lowNum', type=click.INT, default=2, show_default=True, help="max number with low quality allowed")
@click.pass_context
def barcode_pipe(ctx, fq1, fq2, sample, outdir, pattern, whitelist, linkers, lowqual, lownum):
    """
    extract barcode and umi description
    """
    from scopetools.barcode import barcode
    click.echo('barcode pipeline')
    barcode(ctx, fq1, fq2, sample, outdir, pattern, whitelist, linkers, lowqual, lownum)


@cli.command(name='cutadapt', short_help="cutadapt short help")
@click.option('--fq', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True, help="fq help")
@click.option('--sample', type=click.STRING, required=True, help="sample help")
@click.option('--outdir', type=click.Path(file_okay=False, dir_okay=True, writable=True), required=True, help="outdir help")
@click.option('--adapter', cls=MultipleOption, type=AdapterType(), nargs=-1, default=['polyT=A{20}', 'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'], help="adapter help")
@click.option('--minimum-length', type=click.INT, default=20, show_default=True, help="minimum length help")
@click.option('--nextseq-trim', type=click.INT, default=20, show_default=True, help="nextseq trim help")
@click.option('--overlap', type=click.INT, default=5, show_default=True, help="overlap help")
@click.pass_context
def cutadapt_pipe(ctx, fq, sample, outdir, adapter, minimum_length, nextseq_trim, overlap):
    """
    cutadapt description
    """
    from scopetools.cutadapt import cutadapt
    click.echo('cutadapt pipeline')
    cutadapt(ctx, fq, sample, outdir, adapter, minimum_length, nextseq_trim, overlap)


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
    click.echo('STAR pipeline')
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
    click.echo('featureCounts pipeline')
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
    click.echo('count pipeline')
    from scopetools.count import count
    count(ctx, bam, sample, outdir, cells)


@cli.command(name='cluster', short_help="cluster short help")
@click.option('--input', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True, help="scRNA-seq file or path")
@click.option('--outdir', type=click.Path(file_okay=False, dir_okay=True, writable=True), required=True, help="outdir help")
@click.option('--filter-genome', type=click.STRING, default='None', show_default=True, help="filter_genome not analysis")
@click.option('--max-genes', type=click.INT, default=2500, show_default=True, help="the max genes in each cell")
@click.option('--max-permito', type=click.FLOAT, default=0.05, show_default=True, help="the max percent_mito in each cell")
@click.option('--neighbors', type=click.INT, default=15, show_default=True, help="the neighbors")
@click.option('--pc', type=click.INT, default=40, show_default=True, help="Use this many PCs")
@click.option('--rgenes-method', type=click.Choice(['wilcoxon', 't-test', 'logreg']), default='wilcoxon', show_default=True, help="rank genes method,t-test,wilcoxon,logreg")
@click.option('--cluster-algo', type=click.Choice(['leiden', 'louvain']), default='leiden', help="cluster algorithm,leiden or louvain")
@click.option('--plot-method', type=click.Choice(['umap', 'tsne']), default='umap', show_default=True, help="cluster in umap or tsne")
def cluster_pipe(ctx, input, outdir, filter_genome, max_genes, max_permito, neighbors, pc, rgenes_method, cluster_algo, plot_method):
    """
    cluster description
    """
    click.echo('cluster pipeline')
    from scopetools.cluster import cluster
    cluster(ctx, input, outdir, filter_genome, max_genes, max_permito, neighbors, pc, rgenes_method, cluster_algo, plot_method)


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
