# -*- coding: utf-8 -*-
import sys
from pathlib import Path

import pysam

import re
from collections import defaultdict
from .report import Reporter
from .utils import getlogger, CommandWrapper

logger = getlogger(__name__)
logger.setLevel(10)


class GeneName(object):

    def __init__(self, annot):
        self.gtf_file = annot
        self.id_name = defaultdict(str)
        self.parse()

    def parse(self):
        gene_id_pattern = re.compile(r'gene_id "(\S+)";')
        gene_name_pattern = re.compile(r'gene_name "(\S+)"')
        with open(self.gtf_file) as f:
            for line in f.readlines():
                if line.startswith('#!'):
                    continue
                tabs = line.split('\t')
                gtf_type, attributes = tabs[2], tabs[-1]
                if gtf_type == 'gene':
                    gene_id = gene_id_pattern.findall(attributes)[-1]
                    gene_name = gene_name_pattern.findall(attributes)[-1]
                    self.id_name[gene_id] = gene_name

    def __contains__(self, item):
        return item in self.id_name

    def __getitem__(self, item):
        if item in self.id_name:
            return self.id_name[item]
        else:
            return item


class FeatureCountsLogger(object):
    def __init__(self, log, sample):
        self.log = log
        self.stat_info = {
            'visible': {},
            'invisible': {}
        }
        self.parse_log()

    def parse_log(self):
        vals = []
        attrs = ['Assigned', 'Unassigned_NoFeatures', 'Unassigned_Ambiguity']
        with open(self.log, mode='r', encoding='utf-8') as f:
            lines = f.read()

            for p in attrs:
                pattern = re.compile(f'(?<={p})\D*(\d*)')
                vals.append(int(pattern.findall(lines)[0]))
            for attr, val in zip(attrs, vals):
                self.stat_info['visible'][attr] = f'{val} ({val / sum(vals):.2%})'


class Alignment(object):

    def __init__(self, alignment: pysam.AlignedRead):
        self.alignment = alignment
        self.__gene_name = None

    @property
    def gene_id(self):
        if self.alignment.has_tag('XT'):
            return self.alignment.get_tag('XT')
        else:
            return None

    @property
    def gene_name(self):
        return self.__gene_name

    @gene_name.setter
    def gene_name(self, value):
        self.__gene_name = value
        self.alignment.set_tag('XT', self.__gene_name)


def featurecounts(ctx, input, annot, sample, outdir, format, nthreads, type, debug):
    sample_outdir = outdir / sample / '04.featureCounts'
    sample_outdir.mkdir(parents=True, exist_ok=True)

    # featureCounts
    logger.info('featureCounts start!')
    featurecounts_cmd = f'featureCounts -a {annot} -o {sample_outdir / sample} -t {type} -R {format} -T {nthreads} {input}'
    featurecounts_process = CommandWrapper(command=featurecounts_cmd, logger=logger)
    if featurecounts_process.returncode:
        logger.warning('featureCounts error!')
        sys.exit(-1)
    else:
        logger.info('featureCounts done!')

    # samtools sort
    samtools_cmd = f'samtools sort -n -@ {nthreads} -o {sample_outdir / sample}_name_sorted.bam {sample_outdir / input.name}.featureCounts.bam'
    samtools_process = CommandWrapper(samtools_cmd, logger=logger)
    if samtools_process.returncode:
        logger.warning('samtools error!')
        sys.exit(-1)
    else:
        logger.info('samtools done!')

    # convert gene id to gene name in BAM
    gene_name_dict = GeneName(annot)
    with pysam.AlignmentFile(f'{sample_outdir / sample}_name_sorted.bam', mode='rb') as f, pysam.AlignmentFile(f'{sample_outdir / sample}_name_sorted.tmp.bam', mode='wb', template=f) as g:
        for i in f:
            alignment = Alignment(i)
            if alignment.gene_id:
                gene_name = gene_name_dict[alignment.gene_id]
                alignment.gene_name = gene_name
            g.write(alignment.alignment)
    Path(f'{sample_outdir / sample}_name_sorted.tmp.bam').rename(f'{sample_outdir / sample}_name_sorted.bam')

    # parse log
    featurecounts_log = FeatureCountsLogger(log=f'{sample_outdir}/{sample}.summary', sample=sample)

    # report
    logger.info('generate report start!')
    Reporter(name='featureCounts', stat_json=featurecounts_log.stat_info, outdir=sample_outdir.parent)
    logger.info('generate report done!')
