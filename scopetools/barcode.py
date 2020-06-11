# -*- coding: utf-8 -*-
import gzip
import sys

import numpy as np
import pandas as pd
import pysam

from .protocol import BarcodePattern, MisSeq
from .report import Reporter
from .utils import getlogger, CommandWrapper

logger = getlogger(__name__)
logger.setLevel(10)


def barcode(
        ctx,
        fq1s=None,
        fq2s=None,
        sample=None,
        outdir=None,
        bctype=None,
        pattern=None,
        lowqual=None,
        lownum=None,
        whitelist=None,
        linker=None
):
    """

    :param ctx:
    :param fq1s:
    :param fq2s:
    :param sample:
    :param outdir:
    :param bctype:
    :param pattern:
    :param lowqual:
    :param lownum:
    :param whitelist:
    :param linker:
    :return:
    """
    logger.info('Extract barcode start!')
    barcode_pattern = BarcodePattern(pattern)

    cell_len, umi_len = 0, 0
    for start, end in zip(barcode_pattern['C'].start, barcode_pattern['C'].end):
        cell_len += end - start
    for start, end in zip(barcode_pattern['U'].start, barcode_pattern['U'].end):
        umi_len += end - start
    cell_umi_quality_array = np.zeros((cell_len + umi_len, 42), dtype=np.uint64)
    cell_umi_base_array = np.zeros((cell_len + umi_len, 5), dtype=np.uint64)
    base_dict = dict(zip(['A', 'T', 'C', 'G', 'N'], range(0, 5)))

    whitelist_list = []
    if whitelist:
        with open(whitelist, encoding='utf-8', mode='r')as f:
            for line in f.readlines():
                whitelist_list.append(line.strip())
        cell_dict = MisSeq(whitelist_list)
    else:
        cell_dict = None

    if linker:
        with open(linker, encoding='utf-8', mode='r')as f:
            length = [end - start for start, end in zip(barcode_pattern['L'].start, barcode_pattern['L'].end)]
            linker_list = [[] for i in range(len(length))]
            for line in f.readlines():
                linkers = line.strip()
                for nth, val in enumerate(length):
                    linker_list[nth].append(linkers[:val])
                    linkers = linkers[val:]
        linkers_dict = [MisSeq(linker) for linker in linker_list]
    else:
        linkers_dict = None

    sample_outdir = outdir / sample / '01.barcode'
    sample_outdir.mkdir(parents=True, exist_ok=True)
    clean_fastq = sample_outdir / f'{sample}_2.fq.gz'

    if bctype == 'SCOPEv2':
        from .protocol import SCOPEv2 as Sequence
    elif bctype == 'SCOPEv1':
        from .protocol import SCOPEv1 as Sequence
    else:
        from .protocol import Sequence

    with gzip.open(clean_fastq, mode='wt', encoding='utf-8') as f:
        for fq1, fq2 in zip(fq1s, fq2s):
            with pysam.FastxFile(fq1) as f1, pysam.FastxFile(fq2) as f2:
                for seq1, seq2 in zip(f1, f2):
                    if seq1.name == seq2.name:
                        Sequence.seq_info['total_num'] += 1
                    else:
                        logger.warning(f"{fq1} and {fq2} are not compatible")
                        sys.exit(1)

                    sequence = Sequence(
                        seq1=seq1,
                        seq2=seq2,
                        lownum=lownum,
                        lowqual=lowqual,
                        barcode_pattern=barcode_pattern,
                        linkers_dict=linkers_dict,
                        cell_dict=cell_dict
                    )

                    if sequence.rna_sequence:
                        f.write(f'{sequence.rna_sequence}\n')

                        for position, quality in enumerate(sequence.cell_quality + sequence.umi_quality):
                            cell_umi_quality_array[position, quality] += 1

                        for position, base in enumerate(''.join(sequence.cell) + sequence.umi):
                            cell_umi_base_array[position, base_dict[base]] += 1

    # stat
    cell_q30 = cell_umi_quality_array[:cell_len, 30:].sum() / cell_umi_quality_array[:cell_len].sum()
    umi_q30 = cell_umi_quality_array[cell_len:, 30:].sum() / cell_umi_quality_array[cell_len:].sum()

    stat_info = Sequence.stat_info()
    stat_info['Q30 of Barcodes'] = f'{cell_q30:.2%}'
    stat_info['Q30 of UMIs'] = f'{umi_q30:.2%}'

    # indice
    df = pd.DataFrame(cell_umi_base_array)
    df.index += 1
    df.to_csv(
        sample_outdir / 'cell_umi_base.csv',
        mode='w',
        encoding='utf-8',
        header=['A', 'T', 'C', 'G', 'N']
    )

    df = pd.DataFrame(cell_umi_quality_array)
    df.index += 1
    df.to_csv(
        sample_outdir / 'cell_umi_quality.csv',
        mode='w',
        encoding='utf-8',
        header=range(1, cell_len + umi_len)
    )
    logger.info('Extract barcode done!')

    # fastqc
    logger.info('fastqc start!')
    fastqc_cmd = f'fastqc -o {sample_outdir} {clean_fastq}'
    logger.info(fastqc_cmd)
    fastqc_process = CommandWrapper(command=fastqc_cmd, logger=logger)
    if fastqc_process.returncode:
        logger.warning('fastqc error!')
        sys.exit(-1)
    else:
        logger.info('fastqc done!')

    # report
    logger.info('generate report start!')
    Reporter(name='barcode', stat_json=stat_info, outdir=sample_outdir.parent)
    logger.info('generate report done!')
