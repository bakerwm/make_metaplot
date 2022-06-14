#!/usr/bin/env python

"""
Generate bigWig files using deeptools

## Why this script?
In order to generate normalized forward/reverse bigWig files

dUTP library,

# (forward) -–filterRNAstrand=forward keeps minus-strand reads, -f 16
# (reverse) -–filterRNAstrand=reverse keeps plus-strand reads, -F 16

scaleFactor
normalizeUsing 
extendReads
smoothLength 
...
"""

import os
import argparse
from bam2bw import Bam2bw

from utils import (
    make_bam2bw_config, update_obj, load_yaml, dump_yaml, log
)


def show_help(x):
    msg = '\n'.join([
        '-'*80,
        '# 1. Generate a template config file',
        '$ python make_bw.py -t -c {}'.format(x),
        '# 2. Modify the values in YAML' ,
        '# Attentation to the following fields:',
        '  - bam_list: ',
        '  - out_dir: ',
        '  - out_prefix: ',
        '  - normalizeUsing: ',
        '  - binSize: ',
        '  - blackListFileName: ',
        '  - effectiveGenomeSize: ',        
        '# 3. Run the command again',
        '$ python make_bw.py -c {}'.format(x),
        '-'*80,
    ])
    print(msg)


def get_args():
    example = '\n'.join([
        'Example:',
        '# 1. Generate template config file',
        '$ python make_bw.py -c a.yaml -t',
        '# 2: Run program',
        '# modify the config file `a.yaml` according to your data',
        '$ python make_bw.py -c a.yaml',
    ])
    parser = argparse.ArgumentParser(
        prog='make_bw', description='make_bw', epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--config', default=None, required=True,
        help='configs in .yaml file')
    parser.add_argument('-b', '--bam', dest='bam', nargs='+', default=None,
        help='bam files, default: [None]')
    parser.add_argument('-o', dest='out_dir', required=False, default=None,
        help='directory to save bigWig file')
    parser.add_argument('-t', '--get-template', dest='get_template', action='store_true',
        help='Generate the template arguments')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite the plot files')
    return parser


def main():
    args = {
        'bam': None,
        'config': None,
        'get_template': False,
        'overwrite': False,
    }
    args.update(vars(get_args().parse_args()))
    # make sure config.yaml file
    if not isinstance(args['config'], str):
        raise ValueError('config, expect str, got {}'.format(
            type(args['config']).__name__
        ))
    if args['get_template']:
        if os.path.exists(args['config']) and not args['overwrite']:
            log.info('could write to config, file exists: {}'.format(
                args['config']
            ))
        else:
            dump_yaml(make_bam2bw_config(), args['config'])
        show_help(args['config'])
    else:
        d0 = load_yaml(args['config'])
        d1 = make_bam2bw_config()
        d1.update(d0) # from yaml
        # remove None values, from cmd
        if args['bam'] is None:
            args.pop('bam')
        if args['out_dir'] is None:
            args.pop('out_dir')
        d1.update(args) # from cmd
        Bam2bw(**d1).run()
#         print(d1.[''])


if __name__ == '__main__':
    main()

# EOF
