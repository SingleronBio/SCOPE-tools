# -*- coding: utf-8 -*-
from .report import Reporter
from .utils import getlogger

logger = getlogger(__name__)
logger.setLevel(10)


def sample_description(ctx, sample, transcriptome, description, version, outdir, debug):
    sample_outdir = outdir / sample / '00.sample'
    sample_outdir.mkdir(parents=True, exist_ok=True)

    stat_json = {
        'visible': {},
        'invisible': {}
    }
    attrs = [
        'SampleName',
        'Transcriptome',
        'Description',
        'Software Version'
    ]
    vals = [
        sample,
        transcriptome,
        description,
        version
    ]
    for attr, val in zip(attrs, vals):
        stat_json['visible'][attr] = val

    Reporter(name='Sample', stat_json=stat_json, outdir=sample_outdir.parent, debug=debug)
