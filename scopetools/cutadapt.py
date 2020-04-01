# -*- coding: utf-8 -*-
from pathlib import Path
from scopetools.utils import getlogger, CommandWrapper
import json
import re
from scopetools.report import Reporter
import sys

logger = getlogger(__name__)
logger.setLevel(10)


class CutadaptLogger(object):
    def __init__(self, stdout, sample):
        self.stdout = stdout
        self.sample = sample
        self.stat_info = []
        self.parse()

    def parse(self):
        pattern = re.compile(r'(Total reads processed:.*?Total written.*?)\n', flags=re.S)
        pattern_space = re.compile(r'\s{2,}')
        match = pattern.search(self.stdout)
        if match:
            self.stat_info.append(
                {
                    'attr': 'SampleName',
                    'val': self.sample
                }
            )
            for line in match.group().split('\n'):
                if line:
                    line = re.sub(pattern_space, '  ', line)
                    line = line.replace(',', '')
                    attr, val = line.split('  ')
                    self.stat_info.append(
                        {
                            'attr': attr,
                            'val': val
                        }
                    )
        else:
            logger.warning(f'can not match cutadapt log')


def cutadapt(ctx, fq, sample, outdir, adapter, minimum_length, nextseq_trim, overlap):
    adapter_para = ''
    for i, j in zip(['-a'] * len(adapter), adapter):
        adapter_para += f'{i} {j} '
    sample_outdir = Path(outdir, sample, '02.cutadapt')
    sample_outdir.mkdir(parents=True, exist_ok=True)

    # cutadapt
    clean_fastq = sample_outdir / f'{sample}_clean_2.fq.gz'
    cutadapt_cmd = f'cutadapt {adapter_para}-m {minimum_length} --nextseq-trim={nextseq_trim} --overlap {overlap} -o {clean_fastq} {fq}'
    logger.info('cutadapt start!')
    logger.info(cutadapt_cmd)
    cutadapt_process = CommandWrapper(command=cutadapt_cmd, logger=logger)
    if cutadapt_process.returncode:
        logger.warning('cutadapt error!')
        sys.exit(-1)
    else:
        logger.info('cutadapt done!')

    # parse log
    cutadapt_log = CutadaptLogger(stdout=cutadapt_process.stdout, sample=sample)
    with open(sample_outdir / 'stat.json', mode='w', encoding='utf-8') as f:
        json.dump(cutadapt_log.stat_info, f)

    # report
    logger.info('generate report start!')
    Reporter(name='cutadapt', stat_file=sample_outdir / 'stat.json', outdir=sample_outdir.parent)
    logger.info('generate report done!')
