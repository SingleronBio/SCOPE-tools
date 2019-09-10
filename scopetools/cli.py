# -*- coding: utf-8 -*-
import click
from utils import BarcodeType, AdapterType, MultipleOption
from barcode import barcode
from pathlib import Path

__version__ = '0.1.0'
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(chain=True, context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__, prog_name='scopetools')
def cli():
    """Single Cell Omics Preparation Entity Tools"""
    pass


@cli.command(name='barcode', short_help="extract barcode and umi short help")
@click.option('--fq1', type=click.Path(exists=True, file_okay=True, readable=True), required=True, help="read1 fq file")
@click.option('--fq2', type=click.Path(exists=True, file_okay=True, readable=True), required=True, help="read2 fq file")
@click.option('--sample', type=click.STRING, required=True, help="sample name")
@click.option('--outdir', type=click.Path(dir_okay=True, writable=True), required=True, help="output dir")
@click.option('--pattern', type=BarcodeType(), default='C6L15C6L15C6U6T30', show_default=True, help="read1 pattern, C: cellbarcode, L: linker, U: UMI, T: polyT")
@click.option('--whitelist', type=click.Path(exists=True, file_okay=True, readable=True), default=Path(__file__).resolve().parent.joinpath('extra/whitelist/scope/bclist'), show_default=True, help="cell barcode list")
@click.option('--linker1', type=click.Path(exists=True, file_okay=True, readable=True), default=Path(__file__).resolve().parent.joinpath('extra/whitelist/scope/linker1'), show_default=True, help="first linker")
@click.option('--linker2', type=click.Path(exists=True, file_okay=True, readable=True), default=Path(__file__).resolve().parent.joinpath('extra/whitelist/scope/linker2'), show_default=True, help="second linker")
@click.option('--lowQual', type=click.IntRange(0, 30), default=14, show_default=True, help="max phred of base as low quality")
@click.option('--lowNum', type=click.INT, default=2, show_default=True, help="max number with low quality allowed")
@click.pass_context
def barcode_pipe(ctx, fq1, fq2, sample, outdir, pattern, whitelist, linker1, linker2, lowqual, lownum):
    """extract barcode and umi description"""
    click.echo('barcode pipeline')


@cli.command(name='cutadapt', short_help="cutadapt short help")
@click.option('--fq', type=click.Path(exists=True, file_okay=True, readable=True), required=True, help="fq help")
@click.option('--sample', type=click.STRING, required=True, help="sample help")
@click.option('--outdir', type=click.Path(dir_okay=True, writable=True), required=True, help="outdir help")
@click.option('--adapter', cls=MultipleOption, type=AdapterType(), nargs=-1, default=['polyT=A{18}', 'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'], help="adapter help")
@click.option('--minimum-length', type=click.INT, default=20, show_default=True, help="minimum length help")
@click.option('--nextseq-trim', type=click.INT, default=20, show_default=True, help="nextseq trim help")
@click.option('--overlap', type=click.INT, default=5, show_default=True, help="overlap help")
@click.pass_context
def cutadapt_pipe(ctx, fq, sample, outdir, adapter, minimum_length, nextseq_trim, overlap):
    """cutadapt description"""
    click.echo('cutadapt pipeline')


@cli.command(name="STAR", short_help="STAR short help")
@click.option('--fq', type=click.Path(exists=True, file_okay=True, readable=True), required=True, help="fq help")
@click.option('--refFlat', type=click.Path(exists=True, file_okay=True, readable=True), required=True, help="refFlat help")
@click.option('--genomeDir', type=click.Path(dir_okay=True, readable=True), required=True, help="genome help")
@click.option('--sample', type=click.STRING, required=True, help="sample help")
@click.option('--outdir', type=click.Path(dir_okay=True, writable=True), required=True, help="outdir help")
@click.option('--readFilesCommand', type=click.Choice(['zcat', 'bzcat', '-']), default='zcat', show_default=True, help='readFilesCommand help')
@click.option('--runThreadN', type=click.INT, default=2, show_default=True, help='runThreadN help')
@click.pass_context
def star_pipe(ctx, fq, refflat, genomedir, sample, outdir, readfilescommand, runthreadn):
    """STAR description"""
    click.echo('STAR pipeline')


@cli.command(name='featureCounts', short_help="featureCounts short help")
@click.option('--input', type=click.Path(exists=True, file_okay=True, readable=True), required=True, help="input help")
@click.option('--annot', type=click.Path(exists=True, file_okay=True, readable=True), required=True, help="annot help")
@click.option('--sample', type=click.STRING, required=True, help="sample help")
@click.option('--outdir', type=click.Path(dir_okay=True, writable=True), required=True, help="outdir help")
@click.option('--format', type=click.Choice(['SAM', 'BAM', 'CORE']), default='BAM', show_default=True, help="format help")
@click.option('--nthreads', type=click.IntRange(1, 32), default=2, show_default=True, help="nthreads help")
@click.pass_context
def featurecounts_pipe(ctx, input, annot, sample, outdir, format, nthreads):
    """featureCounts description"""
    click.echo('featureCounts pipeline')


@cli.command(name='count', short_help="count short help")
@click.option('--bam', type=click.Path(exists=True, file_okay=True, readable=True), required=True, help="bam help")
@click.option('--sample', type=click.STRING, required=True, help="sample help")
@click.option('--outdir', type=click.Path(dir_okay=True, writable=True), required=True, help="outdir help")
@click.option('--cells', type=click.INT, default=3000, show_default=True, help="sample help")
@click.pass_context
def count_pipe(ctx, bam, sample, outdir, cells):
    """count description"""
    click.echo('count pipeline')


@cli.command(name='run', help="run short help")
@click.option('--fq1', type=click.Path(exists=True, file_okay=True, readable=True), required=True, help="[barcode], read1 fq file")
@click.option('--fq2', type=click.Path(exists=True, file_okay=True, readable=True), required=True, help="[barcode], read2 fq file")
@click.option('--sample', type=click.STRING, required=True, help="[barcode], sample name")
@click.option('--outdir', type=click.Path(dir_okay=True, writable=True), required=True, help="[barcode], output dir")
@click.option('--pattern', type=BarcodeType(), default='C6L15C6L15C6U6T30', show_default=True, help="[barcode], read1 pattern, C: cellbarcode, L: linker, U: UMI, T: polyT")
@click.option('--whitelist', type=click.Path(exists=True, file_okay=True, readable=True), default=Path(__file__).resolve().parent.joinpath('extra/whitelist/scope/bclist'), show_default=True, help="[barcode], cell barcode list")
@click.option('--linker1', type=click.Path(exists=True, file_okay=True, readable=True), default=Path(__file__).resolve().parent.joinpath('extra/whitelist/scope/linker1'), show_default=True, help="[barcode], first linker")
@click.option('--linker2', type=click.Path(exists=True, file_okay=True, readable=True), default=Path(__file__).resolve().parent.joinpath('extra/whitelist/scope/linker2'), show_default=True, help="[barcode], second linker")
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
def run_pipe(ctx, fq1, fq2, sample, outdir, pattern, whitelist, linker1, linker2, lowqual, lownum, adapter, minimum_length, nextseq_trim, overlap, refflat, genomedir, readfilescommand, runthreadn, annot, format, nthreads, cells):
    """run description"""
    click.echo('run pipeline')
    ctx.invoke(barcode_pipe)
    ctx.invoke(cutadapt_pipe)
    ctx.invoke(star_pipe)
    ctx.invoke(featurecounts_pipe)
    ctx.invoke(count_pipe)


if __name__ == '__main__':
    cli()
