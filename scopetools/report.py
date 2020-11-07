# -*- coding: utf-8 -*-
from pathlib import Path

from jinja2 import Environment, select_autoescape, FileSystemLoader

import base64
import json

env = Environment(
    loader=FileSystemLoader(Path(__file__).parent / 'templates'),
    autoescape=select_autoescape(['html', 'xml'])
)


class Reporter(object):
    def __init__(self, name, stat_json: dict, outdir: Path, plot=None, img=None, debug=False):
        self.name = name
        self.stat_json = stat_json
        self.outdir = outdir
        self.plot = plot
        self.img = img
        self.debug = debug
        self.get_report()

    def get_report(self):
        template = env.get_template('base.html')
        json_file = self.outdir / '.data.json'
        if not json_file.exists():
            data = {
                'visible': {},
                'invisible': {}
            }
        else:
            with open(str(json_file), encoding='utf-8', mode='r') as f:
                data = json.load(f)

        data['visible'][self.name + '_summary'] = self.stat_json['visible']
        data['invisible'][self.name + '_summary'] = self.stat_json['invisible']

        if self.plot:
            data['visible'][self.name + '_plot'] = self.plot

        if self.img:
            data['visible'][self.name + '_img'] = {}
            for name in self.img:
                with open(self.img[name], 'rb') as f:
                    data['visible'][self.name + '_img'][name] = base64.b64encode(f.read()).decode()

        with open(str(json_file), encoding='utf-8', mode='w') as f:
            json.dump(data, f)

        with open(str(self.outdir / 'report.html'), encoding='utf-8', mode='w') as f:
            if self.debug:
                data['visible'].update(data['invisible'])
                html = template.render(data['visible'])
            else:
                html = template.render(data['visible'])
            f.write(html)
