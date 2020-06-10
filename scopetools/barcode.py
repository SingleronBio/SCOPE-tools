# -*- coding: utf-8 -*-
import gzip
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Tuple

import numpy as np
import pandas as pd
import pysam

from scopetools.report import Reporter
from scopetools.utils import getlogger, CommandWrapper, cached_property

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


SEQ_INFO = SeqInfo()


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

    def __init__(self, seqs):
        self.seqs: list = seqs
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
        for seq in self.seqs:
            self.raw_seq[seq] = (seq, -1, '', '')

    def gen_seq(self):
        for seq in self.raw_seq:
            self.all_seq[seq] = self.raw_seq[seq]
            for mis in ('A', 'T', 'C', 'G', 'N'):
                for pos, base in enumerate(seq):
                    if mis != base:
                        mis_seq = seq[:pos] + mis + seq[pos + 1:]
                        if mis_seq in self.all_seq:
                            logger.warning(f"{mis_seq}: ({seq}, {pos}, {base}, {mis}) is duplicate with {mis_seq}: {self.all_seq[mis_seq]}")
                        else:
                            self.all_seq[mis_seq] = (seq, pos, base, mis)


class Sequence(object):

    def __init__(
            self,
            seq1: pysam.FastxRecord,
            seq2: pysam.FastxRecord,
            barcode_pattern: BarcodePattern,
            lowqual: int,
            lownum: int
    ):
        self.seq1 = seq1
        self.seq2 = seq2
        self.barcode_pattern = barcode_pattern
        self.lowqual = lowqual
        self.lownum = lownum
        self.stat_info = {}

    @cached_property
    def polyt(self):
        polyt = []
        for start, end in zip(self.barcode_pattern['T'].start, self.barcode_pattern['T'].end):
            polyt.append(self.seq1.sequence[start:end])
        return polyt

    @cached_property
    def cell(self):
        cell = []
        for start, end in zip(self.barcode_pattern['C'].start, self.barcode_pattern['C'].end):
            cell.append(self.seq1.sequence[start:end])
        return cell

    @cached_property
    def cell_quality(self):
        cell_quality = []
        for start, end in zip(self.barcode_pattern['C'].start, self.barcode_pattern['C'].end):
            cell_quality[-1:-1] = [ord(q) - 33 for q in self.seq1.quality[start:end]]
        return cell_quality

    @cached_property
    def umi(self):
        umi = []
        for start, end in zip(self.barcode_pattern['U'].start, self.barcode_pattern['U'].end):
            umi.append(self.seq1.sequence[start:end])
        return ''.join(umi)

    @cached_property
    def umi_quality(self):
        umi_quality = []
        for start, end in zip(self.barcode_pattern['U'].start, self.barcode_pattern['U'].end):
            umi_quality[-1:-1] = [ord(q) - 33 for q in self.seq1.quality[start:end]]
        return umi_quality

    def is_low_quality(self):
        return True if sum(q < self.lowqual for q in self.cell_quality + self.umi_quality) > self.lownum else False

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

    def correct_seq(self):
        pass


class SCOPEv1(Sequence):

    def __init__(self, *args, **kwargs):
        super(SCOPEv1, self).__init__(*args, **kwargs)
        self.rna_sequence = self.correct_seq()

    @cached_property
    def corrected_cell(self):
        return ''.join(self.cell)

    def is_no_polyt(self, strict_t: int = 0, min_t: int = 25):
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

    def is_no_cell(self):
        return ''.join(self.cell).startswith('N')

    def correct_seq(self):
        if self.is_no_polyt():
            SEQ_INFO.no_polyt_num += 1
            return None

        if self.is_low_quality():
            SEQ_INFO.lowqual_num += 1
            return None

        if self.is_no_cell():
            SEQ_INFO.no_cell_num += 1
            return None

        SEQ_INFO.cell_dict[self.corrected_cell] += 1

        # new readID: @barcode_umi_old readID
        rna_sequence = pysam.pysam.FastxRecord()
        rna_sequence.name = f'{self.corrected_cell}_{self.umi}_{self.seq2.name}'
        rna_sequence.sequence = self.seq2.sequence
        rna_sequence.quality = self.seq2.quality
        SEQ_INFO.clean_num += 1
        return rna_sequence


class SCOPEv2(Sequence):

    def __init__(
            self,
            linkers_dict: List[MisSeq],
            cell_dict: MisSeq,
            *args,
            **kwargs
    ):
        super(SCOPEv2, self).__init__(*args, **kwargs)
        self.linkers_dict = linkers_dict
        self.cell_dict = cell_dict
        self._cell = ''
        self.rna_sequence = self.correct_seq()

    @cached_property
    def linkers(self):
        linkers = []
        for start, end in zip(self.barcode_pattern['L'].start, self.barcode_pattern['L'].end):
            linkers.append(self.seq1.sequence[start:end])
        return linkers

    @cached_property
    def corrected_cell(self):
        if 'L' in self.barcode_pattern:
            return self._cell
        else:
            return ''.join(self.cell)

    @cached_property
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

    def correct_seq(self):
        if self.is_no_polyt():
            SEQ_INFO.no_polyt_num += 1
            return None

        if self.is_low_quality():
            SEQ_INFO.lowqual_num += 1
            return None

        if self.is_no_linker():
            SEQ_INFO.no_linker_num += 1
            return None

        if self.is_no_cell():
            SEQ_INFO.no_cell_num += 1
            return None

        if self.corrected_num > 0:
            SEQ_INFO.cell_corrected_num += 1

        SEQ_INFO.cell_dict[self.corrected_cell] += 1

        # new readID: @barcode_umi_old readID
        rna_sequence = pysam.pysam.FastxRecord()
        rna_sequence.name = f'{self.corrected_cell}_{self.umi}_{self.seq2.name}'
        rna_sequence.sequence = self.seq2.sequence
        rna_sequence.quality = self.seq2.quality
        SEQ_INFO.clean_num += 1
        return rna_sequence


def barcode(
        ctx,
        fq1s=None,
        fq2s=None,
        sample=None,
        outdir=None,
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
    :param pattern:
    :param whitelist:
    :param linker:
    :param lowqual:
    :param lownum:
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
    with open(whitelist, encoding='utf-8', mode='r')as f:
        for line in f.readlines():
            whitelist_list.append(line.strip())
    cell_dict = MisSeq(whitelist_list)

    with open(linker, encoding='utf-8', mode='r')as f:
        length = [end - start for start, end in zip(barcode_pattern['L'].start, barcode_pattern['L'].end)]
        linker_list = [[] for i in range(len(length))]
        for line in f.readlines():
            linkers = line.strip()
            for nth, val in enumerate(length):
                linker_list[nth].append(linkers[:val])
                linkers = linkers[val:]
    linkers_dict = [MisSeq(linker) for linker in linker_list]

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

    stat_info = {
        'Number of Reads': SEQ_INFO.total_num,
        'Valid Reads': f'{SEQ_INFO.clean_num} ({SEQ_INFO.clean_num / SEQ_INFO.total_num:.2%})',
        'Valid Barcodes': len(SEQ_INFO.cell_dict),
        'Q30 of Barcodes': f'{cell_q30:.2%}',
        'Q30 of UMIs': f'{umi_q30:.2%}',
        'Reads without polyT': f'{SEQ_INFO.no_polyt_num} ({SEQ_INFO.no_polyt_num / SEQ_INFO.total_num:.2%})',
        'Reads with lowQual': f'{SEQ_INFO.lowqual_num} ({SEQ_INFO.lowqual_num / SEQ_INFO.total_num:.2%})',
        'Reads without linker': f'{SEQ_INFO.no_linker_num} ({SEQ_INFO.no_linker_num / SEQ_INFO.total_num:.2%})',
        'Reads without Barcode': f'{SEQ_INFO.no_cell_num} ({SEQ_INFO.no_cell_num / SEQ_INFO.total_num:.2%})',
        'Reads with corrected Barcode': f'{SEQ_INFO.cell_corrected_num} ({SEQ_INFO.cell_corrected_num / SEQ_INFO.total_num:.2%})',
    }

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
