# -*- coding: utf-8 -*-
import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Tuple

import dnaio

from .utils import cached_property, getlogger

logger = getlogger(__name__)
logger.setLevel(10)


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


@dataclass()
class BarcodeDict(dict):
    start: List[int] = field(default_factory=list)
    end: List[int] = field(default_factory=list)


@dataclass()
class OneSequence(object):
    name: str = ''
    sequence: str = ''
    identifier: str = '+'
    quality: str = ''

    @cached_property
    def rna_sequence(self):
        return f'{self.name}\n{self.sequence}\n{self.identifier}\n{self.quality}'


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
            seq1: dnaio.Sequence,
            seq2: dnaio.Sequence,
            barcode_pattern: BarcodePattern,
            lowqual: int,
            lownum: int,
            linkers_dict: List[MisSeq] = None,
            cell_dict: MisSeq = None
    ):
        self.seq1 = seq1
        self.seq2 = seq2
        self.barcode_pattern = barcode_pattern
        self.lowqual = lowqual
        self.lownum = lownum
        self.linkers_dict = linkers_dict
        self.cell_dict = cell_dict
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
            cell_quality[-1:-1] = [ord(q) - 33 for q in self.seq1.qualities[start:end]]
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
            umi_quality[-1:-1] = [ord(q) - 33 for q in self.seq1.qualities[start:end]]
        return umi_quality


class SCOPEv1(Sequence):
    seq_info = {
        'total_num': 0,
        'clean_num': 0,
        'no_polyt_num': 0,
        'lowqual_num': 0,
        'no_cell_num': 0,
        'cell_dict': defaultdict(int)
    }

    @classmethod
    def add_clean_num(cls, num: int = 1):
        cls.seq_info['clean_num'] += num

    @classmethod
    def add_no_polyt_num(cls, num: int = 1):
        cls.seq_info['no_polyt_num'] += num

    @classmethod
    def add_lowqual_num(cls, num: int = 1):
        cls.seq_info['lowqual_num'] += num

    @classmethod
    def add_no_cell_num(cls, num: int = 1):
        cls.seq_info['no_cell_num'] += num

    @classmethod
    def add_cell_dict(cls, cell, num: int = 1):
        cls.seq_info['cell_dict'][cell] += num

    @classmethod
    def stat_info(cls):
        return {
            'Number of Reads': cls.seq_info['total_num'],
            'Valid Reads': f'{cls.seq_info["clean_num"]} ({cls.seq_info["clean_num"] / cls.seq_info["total_num"]:.2%})',
            'Valid Barcodes': len(cls.seq_info['cell_dict']),
            'Reads without polyT': f'{cls.seq_info["no_polyt_num"]} ({cls.seq_info["no_polyt_num"] / cls.seq_info["total_num"]:.2%})',
            'Reads with lowQual': f'{cls.seq_info["lowqual_num"]} ({cls.seq_info["lowqual_num"] / cls.seq_info["total_num"]:.2%})',
            'Barcode startwith N': f'{cls.seq_info["no_cell_num"]} ({cls.seq_info["no_cell_num"] / cls.seq_info["total_num"]:.2%})',
        }

    def __init__(self, *args, **kwargs):
        super(SCOPEv1, self).__init__(*args, **kwargs)
        self.rna_sequence = self.correct_seq()

    @cached_property
    def corrected_cell(self):
        return ''.join(self.cell)

    @cached_property
    def is_no_polyt(self, strict_t: int = 0, min_t: int = 15):
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

    @cached_property
    def is_no_cell(self):
        return 'N' in ''.join(self.cell)

    @cached_property
    def is_low_quality(self):
        return any(
            [
                sum(q < self.lowqual for q in self.cell_quality) > self.lownum,
                sum(q < self.lowqual for q in self.umi_quality) > self.lownum,
            ]
        )

    def correct_seq(self):
        if self.is_no_polyt:
            self.add_no_polyt_num()
            return None

        if self.is_low_quality:
            self.add_lowqual_num()
            return None

        if self.is_no_cell:
            self.add_no_cell_num()

        self.add_cell_dict(self.corrected_cell)

        # new readID: @barcode_umi_old readID
        rna_sequence = OneSequence()
        rna_sequence.name = f'@{self.corrected_cell}_{self.umi}_{self.seq2.name}'
        rna_sequence.sequence = self.seq2.sequence
        rna_sequence.quality = self.seq2.qualities
        self.add_clean_num()
        return rna_sequence.rna_sequence


class SCOPEv2(Sequence):
    seq_info = {
        'total_num': 0,
        'clean_num': 0,
        'no_polyt_num': 0,
        'lowqual_num': 0,
        'no_linker_num': 0,
        'no_cell_num': 0,
        'cell_corrected_num': 0,
        'cell_dict': defaultdict(int)
    }

    @classmethod
    def add_clean_num(cls, num: int = 1):
        cls.seq_info['clean_num'] += num

    @classmethod
    def add_no_polyt_num(cls, num: int = 1):
        cls.seq_info['no_polyt_num'] += num

    @classmethod
    def add_lowqual_num(cls, num: int = 1):
        cls.seq_info['lowqual_num'] += num

    @classmethod
    def add_no_linker_num(cls, num: int = 1):
        cls.seq_info['no_linker_num'] += num

    @classmethod
    def add_no_cell_num(cls, num: int = 1):
        cls.seq_info['no_cell_num'] += num

    @classmethod
    def add_cell_corrected_num(cls, num: int = 1):
        cls.seq_info['cell_corrected_num'] += num

    @classmethod
    def add_cell_dict(cls, cell, num: int = 1):
        cls.seq_info['cell_dict'][cell] += num

    @classmethod
    def stat_info(cls):
        return {
            'Number of Reads': cls.seq_info['total_num'],
            'Valid Reads': f'{cls.seq_info["clean_num"]} ({cls.seq_info["clean_num"] / cls.seq_info["total_num"]:.2%})',
            'Valid Barcodes': len(cls.seq_info['cell_dict']),
            'Reads without polyT': f'{cls.seq_info["no_polyt_num"]} ({cls.seq_info["no_polyt_num"] / cls.seq_info["total_num"]:.2%})',
            'Reads with lowQual': f'{cls.seq_info["lowqual_num"]} ({cls.seq_info["lowqual_num"] / cls.seq_info["total_num"]:.2%})',
            'Reads without linker': f'{cls.seq_info["no_linker_num"]} ({cls.seq_info["no_linker_num"] / cls.seq_info["total_num"]:.2%})',
            'Reads without Barcode': f'{cls.seq_info["no_cell_num"]} ({cls.seq_info["no_cell_num"] / cls.seq_info["total_num"]:.2%})',
            'Reads with corrected Barcode': f'{cls.seq_info["cell_corrected_num"]} ({cls.seq_info["cell_corrected_num"] / cls.seq_info["total_num"]:.2%})',
        }

    def __init__(
            self,
            *args,
            **kwargs
    ):
        super(SCOPEv2, self).__init__(*args, **kwargs)
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
        return self._cell

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

    @cached_property
    def is_no_linker(self):
        return any(
            [
                linker not in linker_dict for linker, linker_dict in zip(self.linkers, self.linkers_dict)
            ]
        )

    @cached_property
    def is_no_cell(self):
        if self.corrected_num < -10 or self.corrected_num > 1:
            return True
        else:
            return False

    @cached_property
    def is_low_quality(self):
        return True if sum(q < self.lowqual for q in self.cell_quality + self.umi_quality) > self.lownum else False

    @cached_property
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
        if self.is_no_polyt:
            self.add_no_polyt_num()
            return None

        if self.is_low_quality:
            self.add_lowqual_num()
            return None

        if self.is_no_linker:
            self.add_no_linker_num()
            return None

        if self.is_no_cell:
            self.add_no_cell_num()
            return None

        if self.corrected_num > 0:
            self.add_cell_corrected_num()

        self.add_cell_dict(self.corrected_cell)

        # new readID: @barcode_umi_old readID
        rna_sequence = OneSequence()
        rna_sequence.name = f'@{self.corrected_cell}_{self.umi}_{self.seq2.name}'
        rna_sequence.sequence = self.seq2.sequence
        rna_sequence.quality = self.seq2.qualities
        self.add_clean_num()
        return rna_sequence.rna_sequence
