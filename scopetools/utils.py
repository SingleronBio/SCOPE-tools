# -*- coding: utf-8 -*-


import click
import re
import logging
import sys


class BarcodeType(click.ParamType):
    name = "barcode pattern"

    def convert(self, value, param, ctx):
        p = re.compile(r'([CLTU]\d+)+')
        if re.fullmatch(p, value):
            return value
        else:
            click.echo('fail')
            raise click.BadParameter("{} is not a valid adapter pattern".format(value))


class AdapterType(click.ParamType):
    name = "adapter pattern"

    def convert(self, value, param, ctx):
        """
        regular expressions match the following pattern
        polyT=A{5}T{10}
        p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
        A{18}
        AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
        poly=a{18}
        p5=agatcagatcagatc
        a{18}
        agatcagatcagatc
        """
        poly_pattern = re.compile(r'(^[ATCGatcg]{\d+})+$')
        p5_pattern = re.compile(r'(^[ATCGatcg]+)$')
        adapters = []
        for i in value:
            adapter = i.split('=')[-1]
            if re.fullmatch(poly_pattern, adapter) or re.fullmatch(p5_pattern, adapter):
                adapters.append(adapter)
            else:
                click.echo('fail')
                raise click.BadParameter("{} is not a valid adapter pattern".format(adapter))
                # self.fail("{} is not a valid adapter pattern".format(adapter), param, ctx)
        return adapters


class MultipleOption(click.Option):

    def __init__(self, *args, **kwargs):
        self.save_other_options = kwargs.pop('save_other_options', True)
        nargs = kwargs.pop('nargs', -1)
        assert nargs == -1, 'nargs, if set, must be -1 not {}'.format(nargs)
        super(MultipleOption, self).__init__(*args, **kwargs)
        self._previous_parser_process = None
        self._eat_all_parser = None

    def add_to_parser(self, parser, ctx):

        def parser_process(value, state):
            # method to hook to the parser.process
            done = False
            value = [value]
            if self.save_other_options:
                # grab everything up to the next option
                while state.rargs and not done:
                    for prefix in self._eat_all_parser.prefixes:
                        if state.rargs[0].startswith(prefix):
                            done = True
                    if not done:
                        value.append(state.rargs.pop(0))
            else:
                # grab everything remaining
                value += state.rargs
                state.rargs[:] = []
            value = tuple(value)

            # call the actual process
            self._previous_parser_process(value, state)

        retval = super(MultipleOption, self).add_to_parser(parser, ctx)
        for name in self.opts:
            our_parser = parser._long_opt.get(name) or parser._short_opt.get(name)
            if our_parser:
                self._eat_all_parser = our_parser
                self._previous_parser_process = our_parser.process
                our_parser.process = parser_process
                break
        return retval


class MutuallyExclusiveOption(click.Option):
    def __init__(self, *args, **kwargs):
        self.mutually_exclusive = set(kwargs.pop('mutually_exclusive', []))
        help = kwargs.get('help', '')
        if self.mutually_exclusive:
            ex_str = ', '.join(self.mutually_exclusive)
            kwargs['help'] = help + (
                    ' NOTE: This argument is mutually exclusive with '
                    ' arguments: [' + ex_str + '].'
            )
        super(MutuallyExclusiveOption, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        if self.mutually_exclusive.intersection(opts) and self.name in opts:
            raise click.UsageError(
                "Illegal usage: `{}` is mutually exclusive with arguments `{}`.".format(self.name, ', '.join(self.mutually_exclusive))
            )

        return super(MutuallyExclusiveOption, self).handle_parse_result(ctx, opts, args)


class SingleLevelFilter(logging.Filter):
    """
    reject is True, only include passlevel message
    reject is False, exclude passlevel message
    """

    def __init__(self, passlevel: int, reject: bool):
        super().__init__()
        self.passlevel = passlevel
        self.reject = reject

    def filter(self, record):
        if self.reject:
            return record.levelno != self.passlevel
        else:
            return record.levelno == self.passlevel


def getlogger(name):
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    data_fmt = '%Y-%M-%D %H:%M:%S'
    formatter = logging.Formatter(log_fmt, data_fmt)
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    stdout = logging.StreamHandler(sys.stdout)
    stdout.setFormatter(formatter)
    stdout.addFilter(SingleLevelFilter(logging.INFO, False))
    logger.addHandler(stdout)
    stderr = logging.StreamHandler(sys.stderr)
    stderr.setFormatter(formatter)
    stderr.addFilter(SingleLevelFilter(logging.WARN, False))
    logger.addHandler(stderr)
    return logger
