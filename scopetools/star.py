# -*- coding: utf-8 -*-
import re
import sys
from pathlib import Path

from scopetools.report import Reporter
from scopetools.utils import getlogger, CommandWrapper

logger = getlogger(__name__)
logger.setLevel(10)


class STARLogger(object):
    def __init__(self, star_log, picard_log, sample):
        self.star_log = star_log
        self.picard_log = picard_log
        self.stat_info = {
            'SampleName': f'{sample}'
        }
        self.plot = {}
        self.parse_star()
        self.parse_picard()

    def parse_star(self):
        with open(self.star_log, mode='r', encoding='utf-8') as f:
            lines = f.read()
            uniquely_pattern = re.compile(r'Uniquely mapped reads\D*(\d*\.?\d*%?)', flags=re.S)
            multiple_pattern = re.compile(r'of reads mapped to too many loci\D*(\d*\.?\d*%?)', flags=re.S)
            i, j = uniquely_pattern.findall(lines)
            self.stat_info['Uniquely_mapped'] = f'{i} ({j})'
            i, j = multiple_pattern.findall(lines)
            self.stat_info['Multiple_mapped'] = f'{i} ({j})'

    def parse_picard(self):
        with open(self.picard_log, mode='r', encoding='utf-8') as f:
            lines = f.read()
            metrics_pattern = re.compile(r'(?<=picard.analysis.RnaSeqMetrics)\n(.*)\n(.*)\n')
            match = metrics_pattern.search(lines)
            if match:
                attrs = match.groups()[0].split('\t')
                vals = match.groups()[1].split('\t')
                picard_dict = dict(zip(attrs, vals))
                self.stat_info['Exonic_Regions'] = f'{int(picard_dict["CODING_BASES"]) + int(picard_dict["UTR_BASES"]):d} ({float(picard_dict["PCT_CODING_BASES"]) + float(picard_dict["PCT_UTR_BASES"]):.2%})'
                self.stat_info['Intronic_Regions'] = f'{int(picard_dict["INTRONIC_BASES"]):d} ({float(picard_dict["PCT_INTRONIC_BASES"]):.2%})'
                self.stat_info['Intergenic_Regions'] = f'{int(picard_dict["INTERGENIC_BASES"]):d} ({float(picard_dict["PCT_INTRONIC_BASES"]):.2%})'
                self.plot = {
                    'region_labels': ['Exonic Regions', 'Intronic Regions', 'Intergenic Regions'],
                    'region_values': [int(picard_dict["CODING_BASES"]), int(picard_dict["INTRONIC_BASES"]), int(picard_dict["INTERGENIC_BASES"])]
                }


def star(ctx, fq, refflat, genomedir, sample, outdir, readfilescommand, runthreadn):
    sample_outdir = Path(outdir, sample, '03.STAR')
    sample_outdir.mkdir(parents=True, exist_ok=True)

    # STAR
    out_prefix = sample_outdir / f'{sample}_'
    star_cmd = f'STAR --runThreadN {runthreadn} --genomeDir {genomedir} --readFilesIn {fq} --readFilesCommand {readfilescommand} --outFilterMultimapNmax 1 --outFileNamePrefix {out_prefix} --outSAMtype BAM SortedByCoordinate'
    logger.info('STAR start!')
    star_process = CommandWrapper(star_cmd, logger=logger)
    if star_process.returncode:
        logger.warning('STAR error!')
        sys.exit(-1)
    else:
        logger.info('STAR done')

    # picard
    out_bam = sample_outdir / f'{sample}_Aligned.sortedByCoord.out.bam'
    region_txt = sample_outdir / f'{sample}_region.log'
    picard_cmd = f'picard -Xmx4G -XX:ParallelGCThreads=4 CollectRnaSeqMetrics I={out_bam} O={region_txt} REF_FLAT={refflat} STRAND=NONE VALIDATION_STRINGENCY=SILENT'
    logger.info('stat mapping region start!')
    picard_process = CommandWrapper(picard_cmd, logger=logger)
    if picard_process.returncode:
        logger.warning('stat mapping region error!')
        sys.exit(-1)
    else:
        logger.info('stat mapping region done!')

    # parse_log
    star_log = STARLogger(star_log=f'{out_prefix}Log.final.out', picard_log=region_txt, sample=sample)

    # report
    logger.info('generate report start!')
    Reporter(name='STAR', stat_json=star_log.stat_info, outdir=sample_outdir.parent, plot=star_log.plot)

    logger.info('generate report done!')
