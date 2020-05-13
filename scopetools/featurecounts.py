# -*- coding: utf-8 -*-
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import pysam

from scopetools.report import Reporter
from scopetools.utils import getlogger, CommandWrapper

logger = getlogger(__name__)
logger.setLevel(10)


class FeatureCountsLogger(object):
    def __init__(self, log, sample):
        self.log = log
        self.stat_info = [
            {
                'attr': 'SampleName',
                'val': f'{sample}'
            }
        ]
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
                self.stat_info.append(
                    {
                        'attr': attr,
                        'val': f'{val} ({val / sum(vals):.2%})'
                    }
                )


def featurecounts(ctx, input, annot, sample, outdir, format, nthreads):
    sample_outdir = Path(outdir, sample, '04.featureCounts')
    sample_outdir.mkdir(parents=True, exist_ok=True)

    # featureCounts
    logger.info('featureCounts start!')
    featurecounts_cmd = f'featureCounts -a {annot} -o {sample_outdir / sample} -R {format} -T {nthreads} {input}'
    featurecounts_process = CommandWrapper(command=featurecounts_cmd, logger=logger)
    if featurecounts_process.returncode:
        logger.warning('featureCounts error!')
        sys.exit(-1)
    else:
        logger.info('featureCounts done!')

        # samtools sort
    samtools_cmd = f'samtools sort -n -@ {nthreads} -o {sample_outdir / sample}_name_sorted.bam {sample_outdir / Path(input).name}.featureCounts.bam'
    samtools_process = CommandWrapper(samtools_cmd, logger=logger)
    if samtools_process.returncode:
        logger.warning('samtools error!')
        sys.exit(-1)
    else:
        logger.info('samtools done!')

    # parse log
    featurecounts_log = FeatureCountsLogger(log=f'{sample_outdir}/{sample}.summary', sample=sample)
    with open(sample_outdir / 'stat.json', mode='w', encoding='utf-8') as f:
        json.dump(featurecounts_log.stat_info, f)

    # report
    logger.info('generate report start!')
    Reporter(name='featureCounts', stat_file=sample_outdir / 'stat.json', outdir=sample_outdir.parent)
    logger.info('generate report done!')
