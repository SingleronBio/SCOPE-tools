# -*- coding: utf-8 -*-
import re
from collections import defaultdict
import numpy as np
import pandas as pd
import pysam
import sys
from pathlib import Path
import gzip
from scopetools.utils import getlogger, CommandWrapper
from dataclasses import dataclass, field
from typing import List, Dict, Tuple
import json
from scopetools.report import Reporter

logger = getlogger(__name__)
logger.setLevel(10)


@dataclass()
class SeqInfo(object):
    total_num: int = 0
    clean_num: int = 0
    no_polyt_num: int = 0
    lowqual_num: int = 0
    no_linker_num: int = 0
    no_cell_num: int = 0
    cell_corrected_num: int = 0
    cell_dict = defaultdict(int)


@dataclass()
class BarcodeDict(dict):
    start: List[int] = field(default_factory=list)
    end: List[int] = field(default_factory=list)


class BarcodePattern(object):

    def __init__(self, pattern):
        self.pattern = pattern
        self.barcode_dict = defaultdict(BarcodeDict)
        self.parse()

    def parse(self):
        start = 0
        p = re.findall(r'([CLNTU])(\d+)', self.pattern)
        for i in p:
            end = start + int(i[1])
            self.barcode_dict[i[0]].start.append(start)
            self.barcode_dict[i[0]].end.append(end)
            start = end

    def __getitem__(self, key):
        return self.barcode_dict[key]

    def __contains__(self, key):
        if key in self.barcode_dict:
            return True
        else:
            return False


class MisSeq(object):

    def __init__(self, file):
        self.file = file
        self.raw_seq: Dict[str, Tuple[str, int, str, str]] = {}
        self.all_seq: Dict[str, Tuple[str, int, str, str]] = {}
        self.read_seq()
        self.gen_seq()

    def __getitem__(self, seq: str):
        if seq in self.all_seq:
            return self.all_seq[seq]
        else:
            return None

    def __contains__(self, seq: str):
        return seq in self.all_seq

    def read_seq(self):
        with open(self.file, mode='r', encoding='utf-8') as f:
            for line in f.readlines():
                seq = line.strip()
                if seq:
                    self.raw_seq[seq] = (seq, -1, '', '')

    def gen_seq(self):
        for seq in self.raw_seq:
            self.all_seq[seq] = self.raw_seq[seq]
            for mis in ('A', 'T', 'C', 'G'):
                for pos, base in enumerate(seq):
                    if mis != base:
                        mis_seq = seq[:pos] + mis + seq[pos + 1:]
                        if mis_seq in self.all_seq:
                            logger.warning(f"{mis_seq}: ({seq}, {pos}, {base}, {mis}) is duplicate with {mis_seq}: {self.all_seq[mis_seq]}")
                        else:
                            self.all_seq[mis_seq] = (seq, pos, base, mis)


class Sequence(object):

    def __init__(self,
                 seq1: pysam.FastxRecord,
                 seq2: pysam.FastxRecord,
                 barcode_pattern: BarcodePattern,
                 linkers_dict: List[MisSeq],
                 cell_dict: MisSeq):
        """

        :param seq1:
        :param seq2:
        :param barcode_pattern:
        :param linkers_dict:
        :param cell_dict:
        """
        self.seq1 = seq1
        self.seq2 = seq2
        self.barcode_pattern = barcode_pattern
        self.linkers_dict = linkers_dict
        self.cell_dict = cell_dict
        self._cell = ''

    @property
    def linkers(self):
        linkers = []
        for start, end in zip(self.barcode_pattern['L'].start, self.barcode_pattern['L'].end):
            linkers.append(self.seq1.sequence[start:end])
        return linkers

    @property
    def polyt(self):
        polyt = []
        for start, end in zip(self.barcode_pattern['T'].start, self.barcode_pattern['T'].end):
            polyt.append(self.seq1.sequence[start:end])
        return polyt

    @property
    def cell(self):
        cell = []
        for start, end in zip(self.barcode_pattern['C'].start, self.barcode_pattern['C'].end):
            cell.append(self.seq1.sequence[start:end])
        return cell

    @property
    def cell_quality(self):
        cell_quality = []
        for start, end in zip(self.barcode_pattern['C'].start, self.barcode_pattern['C'].end):
            cell_quality[-1:-1] = [ord(q) - 33 for q in self.seq1.quality[start:end]]
        return cell_quality

    @property
    def umi(self):
        umi = []
        for start, end in zip(self.barcode_pattern['U'].start, self.barcode_pattern['U'].end):
            umi.append(self.seq1.sequence[start:end])
        return ''.join(umi)

    @property
    def umi_quality(self):
        umi_quality = []
        for start, end in zip(self.barcode_pattern['U'].start, self.barcode_pattern['U'].end):
            umi_quality[-1:-1] = [ord(q) - 33 for q in self.seq1.quality[start:end]]
        return umi_quality

    @property
    def corrected_cell(self):
        if 'L' in self.barcode_pattern:
            return self._cell
        else:
            return ''.join(self.cell)

    @property
    def corrected_num(self):
        corrected_num = 0
        for cell in self.cell:
            if cell in self.cell_dict:
                self._cell += self.cell_dict[cell][0]
                if self.cell_dict[cell][1] == -1:
                    corrected_num += 0
                else:
                    corrected_num += 1
                    logger.warning(f"cell barcode {cell} is corrected as {self.cell_dict[cell][0]}")
            else:
                corrected_num += -100
        return corrected_num

    def is_low_quality(self, lowqual: int, lownum: int):
        """

        :param lowqual:
        :param lownum:
        :return:
        """
        return True if sum(q < lowqual for q in self.cell_quality + self.umi_quality) > lownum else False

    def is_no_polyt(self, strict_t: int = 0, min_t: int = 10):
        """

        :param strict_t:
        :param min_t:
        :return:
        """
        polyt = ''.join(self.polyt)
        return any(
            [
                polyt[:strict_t].count('T') < strict_t,
                polyt.count('T') < min_t
            ]
        )

    def is_no_linker(self):
        return any(
            [
                linker not in linker_dict for linker, linker_dict in zip(self.linkers, self.linkers_dict)
            ]
        )

    def is_no_cell(self):
        if self.corrected_num < -10 or self.corrected_num > 1:
            return True
        else:
            return False


def barcode(ctx, fq1, fq2, sample, outdir, pattern, whitelist, linkers, lowqual, lownum):
    """

    :param ctx:
    :param fq1:
    :param fq2:
    :param sample:
    :param outdir:
    :param pattern:
    :param whitelist:
    :param linkers:
    :param lowqual:
    :param lownum:
    :return:
    """
    logger.info('Extract barcode start!')
    barcode_pattern = BarcodePattern(pattern)
    seq_info = SeqInfo()

    cell_len, umi_len = 0, 0
    for start, end in zip(barcode_pattern['C'].start, barcode_pattern['C'].end):
        cell_len += end - start
    for start, end in zip(barcode_pattern['U'].start, barcode_pattern['U'].end):
        umi_len += end - start
    cell_umi_quality_array = np.zeros((cell_len + umi_len, 42), dtype=np.uint64)
    cell_umi_base_array = np.zeros((cell_len + umi_len, 5), dtype=np.uint64)
    base_dict = dict(zip(['A', 'T', 'C', 'G', 'N'], range(0, 5)))

    linkers_dict = [MisSeq(linker) for linker in linkers]
    cell_dict = MisSeq(whitelist)

    sample_outdir = Path(outdir, sample, '01.barcode')
    sample_outdir.mkdir(parents=True, exist_ok=True)
    clean_fastq = sample_outdir / f'{sample}_2.fq.gz'
    with pysam.FastxFile(fq1) as f1, pysam.FastxFile(fq2) as f2, gzip.open(clean_fastq, mode='wt', encoding='utf-8') as f:
        for seq1, seq2 in zip(f1, f2):
            if seq1.name == seq2.name:
                seq_info.total_num += 1
            else:
                logger.warning(f"{fq1} and {fq2} are not compatible")
                sys.exit(1)

            sequence = Sequence(seq1, seq2, barcode_pattern, linkers_dict, cell_dict)
            if 'T' in barcode_pattern:
                if sequence.is_no_polyt():
                    seq_info.no_polyt_num += 1
                    continue

            if sequence.is_low_quality(lowqual=lowqual, lownum=lownum):
                seq_info.lowqual_num += 1
                continue

            if sequence.is_no_linker():
                seq_info.no_linker_num += 1
                continue

            if sequence.is_no_cell():
                seq_info.no_cell_num += 1
                continue

            if sequence.corrected_num > 0:
                seq_info.cell_corrected_num += 1

            seq_info.cell_dict[sequence.corrected_cell] += 1

            # new readID: @barcode_umi_old readID
            rna_sequence = pysam.FastxRecord()
            rna_sequence.name = f'{sequence.corrected_cell}_{sequence.umi}_{seq2.name}'
            rna_sequence.sequence = seq2.sequence
            rna_sequence.quality = seq2.quality
            f.write(f'{rna_sequence}\n')
            seq_info.clean_num += 1

            for position, quality in enumerate(sequence.cell_quality + sequence.umi_quality):
                cell_umi_quality_array[position, quality] += 1

            for position, base in enumerate(''.join(sequence.cell) + sequence.umi):
                cell_umi_base_array[position, base_dict[base]] += 1

    # stat
    cell_q30 = cell_umi_quality_array[:cell_len, 30:].sum() / cell_umi_quality_array[:cell_len].sum()
    umi_q30 = cell_umi_quality_array[cell_len:, 30:].sum() / cell_umi_quality_array[cell_len:].sum()

    with open(sample_outdir / 'stat.json', mode='w', encoding='utf-8') as f:
        stat_info = [
            {
                'attr': 'SampleName',
                'val': sample
            },
            {
                'attr': 'Number of Reads',
                'val': seq_info.total_num,
            },
            {
                'attr': 'Valid Reads',
                'val': f'{seq_info.clean_num} ({seq_info.clean_num / seq_info.total_num:.2%})',
            },
            {
                'attr': 'Valid Barcodes',
                'val': len(seq_info.cell_dict),
            },
            {
                'attr': 'Q30 of Barcodes',
                'val': f'{cell_q30:.2%}',
            },
            {
                'attr': 'Q30 of UMIs',
                'val': f'{umi_q30:.2%}',
            },
            {
                'attr': 'Reads without polyT',
                'val': f'{seq_info.no_polyt_num} ({seq_info.no_polyt_num / seq_info.total_num:.2%})',
            },
            {
                'attr': 'Reads with lowQual',
                'val': f'{seq_info.lowqual_num} ({seq_info.lowqual_num / seq_info.total_num:.2%})',
            },
            {
                'attr': 'Reads without linker',
                'val': f'{seq_info.no_linker_num} ({seq_info.no_linker_num / seq_info.total_num:.2%})',
            },
            {
                'attr': 'Reads without Barcode',
                'val': f'{seq_info.no_cell_num} ({seq_info.no_cell_num / seq_info.total_num:.2%})',
            },
            {
                'attr': 'Reads with corrected Barcode',
                'val': f'{seq_info.cell_corrected_num} ({seq_info.cell_corrected_num / seq_info.total_num:.2%})'
            }
        ]
        json.dump(stat_info, f)

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
        sample_outdir / 'cell_umi_qualitye.csv',
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
    Reporter(name='barcode', stat_file=sample_outdir / 'stat.json', outdir=sample_outdir.parent)
    logger.info('generate report done!')
