# -*- coding: utf-8 -*-
from pathlib import Path
import json
import pandas as pd
from jinja2 import Environment, PackageLoader, select_autoescape, FileSystemLoader

env = Environment(
    loader=FileSystemLoader(Path(__file__).parent / 'templates'),
    autoescape=select_autoescape(['html', 'xml'])
)


class Reporter(object):
    def __init__(self, name, stat_file, outdir: Path, plot=None):
        self.name = name
        self.stat_file = stat_file
        self.outdir = outdir
        self.plot = plot
        self.get_report()

    def get_report(self):
        template = env.get_template('base.html')
        json_file = self.outdir / '.data.json'
        if not json_file.exists():
            data = {}
        else:
            with open(json_file, encoding='utf-8', mode='r') as f:
                data = json.load(f)

        df = pd.read_json(self.stat_file)
        data[self.name + '_summary'] = df.values.tolist()

        if self.plot:
            data[self.name + '_plot'] = self.plot

        with open(self.outdir / 'report.html', encoding='utf-8', mode='w') as f:
            html = template.render(data)
            f.write(html)

        with open(json_file, encoding='utf-8', mode='w') as f:
            json.dump(data, f)
