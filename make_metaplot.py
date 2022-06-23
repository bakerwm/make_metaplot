#!/usr/bin/env python
#-*- encoding:utf-8 -*-
"""
Generate metaplot using deeptools toolkit

## Why this script?
In order to generate strand-specific genomebody coverage plots for RNAseq data

    Sense: bigWig(fwd) gene(+); bigWig(rev) gene(-)
Antisense: bigWig(fwd) gene(-); bigWig(fwd) gene(+)

## How-To

### A. DNA seq data (non-stranded)
# eg: ChIP-Seq, ATAC-seq, CUT&Tag, ...
1. BAM -> bigWig
2. bigWig + BED -> matrix
3. matrix -> plots (metaplot, heatmap, ...)

### B. RNA seq data (stranded)
# eg: RNA-Seq, ChrRNAseq, TTseq, ...
1. BAM -> bigWig (fwd, rev)
2. bigWig(fwd) + BED(+) -> matrix-1 (fwd)
3. bigWig(rev) + BED(-) -> matrix-2 (fwd)
4. bigWig(fwd) + BED(-) -> matrix-3 (rev)
5. bigWig(rev) + BED(+) -> matrix-4 (rev)
6. merge matrix (by rows)
   matrix-1 + matrix-2 -> matrix_sens
   matrix-3 + matrix-4 -> matrix_anti
7. matrix (sens/anti) -> plots
"""

import os
import argparse
from bam2profile import Bam2profile
from bam2heatmap import Bam2heatmap
from bw2profile import Bw2profile
from bw2heatmap import Bw2heatmap
from matrix2profile import Matrix2profile
from matrix2heatmap import Matrix2heatmap
from parse_args import get_init_args
from utils import (
    load_yaml, dump_yaml, is_valid_bam, is_valid_bigwig, is_valid_matrix, log,
    is_valid_file, get_program_msg
)


def func_to_plot(plot_type='profile', **kwargs):
    """
    priority:
    1. Matrix2profile: matrix
    2. Bw2profile: bw_fwd_list, bw_rev_list
    3. Bw2profile: bw_list
    4. Bam2profile: bam_list
    5. Matrix2profile: matrix
    """
    args = kwargs
    # variables
    m = args.get('matrix', None)
    bf = args.get('bw_fwd_list', None)
    br = args.get('bw_rev_list', None)
    bw = args.get('bw_list', None)
    bam = args.get('bam_list', None)
    # 1. matrix
    if is_valid_matrix(m):
        fun1 = Matrix2profile
        fun2 = Matrix2heatmap
        args.update({
            'bw_fwd_list': None,
            'bw_rev_list': None,
            'bw_list': None,
            'bam_list': None,
        })
    # elif is_valid_file(bf, is_valid_bigwig) and is_valid_file(br, is_valid_bigwig):
    elif is_valid_file([bf, br], is_valid_bigwig):
        fun1 = Bw2profile
        fun2 = Bw2heatmap
        args.update({
            'matrix': None,
            'bw_list': None,
            'bam_list': None,
        })
    elif is_valid_file(bw, is_valid_bigwig):
        fun1 = Bw2profile
        fun2 = Bw2heatmap
        args.update({
            'matrix': None,
            'bw_fwd_list': None,
            'bw_rev_list': None,
            'bam_list': None,
        })
    elif is_valid_file(bam, is_valid_bam):
        fun1 = Bam2profile
        fun2 = Bam2heatmap
        args.update({
            'matrix': None,
            'bw_fwd_list': None,
            'bw_rev_list': None,
            'bw_list': None,
        })
    else:
        fun1 = fun2 = None
    fun = fun1 if plot_type == 'profile' else fun2
    return (fun, args)


def show_help(x):
    msg = '\n'.join([
        '='*80,
        '# 1. Generate a template config file',
        '$ python make_metaplot.py -t -c {}'.format(x),
        '# 2. Modify the values in YAML' ,
        '# Attentation to the following fields:',
        '  - bam_list: ',
        '  - bw_list: ',
        '  - bw_fwd_list: ',
        '  - bw_rev_list: ',
        '  - matrix: ',
        '  - region_list: ',
        '  - samplesLabel: sample1 sample2',
        '  - regionsLabel: gene',
        '  - colors red black',
        '  - colorList: red,white,blue',
        '  - colorMap: Reds',
        '  - yMin, yMax',
        '# 3. Run the command again',
        '$ python make_profile.py -c {}'.format(x),
        '='*80,
    ])
    print(msg)


def get_args():
    example = '\n'.join([
        'Example:',
        '# 1. Generate template config file',
        '$ python make_metaplot.py -c a.yaml -t',
        '# 2: Run program',
        '# modify the config file `a.yaml` according to your data',
        '$ python make_metaplot.py -c a.yaml',
    ])
    parser = argparse.ArgumentParser(
        prog='make_metaplot', description='make_metaplot', epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--config', default=None, required=True,
        help='configs in .yaml file')
    parser.add_argument('-t', '--get-template', dest='get_template', action='store_true',
        help='Generate the template arguments')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite the plot files')
    return parser


def main():
    args = get_args().parse_args() # input
    # check config
    args_init = get_init_args() #
    if args.get_template:
        try:
            dump_yaml(args_init, args.config)
        except IOError as e:
            log.error(e)
        show_help(args.config)
    else:
        d1 = load_yaml(args.config) # i/o config
        if not isinstance(d1, dict):
            log.error('Coult not read config file: {}'.format(args.config))
            return None # skipped
        args_init.update(d1)
        fun1, args1 = func_to_plot('profile', **args_init) # for profile
        fun2, args2 = func_to_plot('heatmap', **args_init) # for heatmap
        if fun1 is None or fun2 is None:
            log.error('Could not determine the program, check config file')
            return None
        # show messages
        msg = get_program_msg(**d1) # input
        print(msg)
        fun1(**args1).run()
        fun2(**args2).run()

    # # make sure config.yaml file
    # if not isinstance(args['config'], str):
    #     raise ValueError('config, expect str, got {}'.format(
    #         type(args['config']).__name__
    #     ))
    # if args['get_template']:
    #     if os.path.exists(args['config']) and not args['overwrite']:
    #         log.info('could write to config, file exists: {}'.format(
    #             args['config']
    #         ))
    #     else:
    #         dump_yaml(make_config(), args['config'])
    #     show_help(args['config'])
    # else:
    #     d1 = make_config()
    #     if d1 is None:
    #         print('error, check config file: {}'.format(args['config']))
    #         sys.exit(1)
    #     d1.update(load_yaml(args['config']))
    #     fun, d2 = find_fun(d1)
    #     if fun is None:
    #         raise ValueError('unable to determin functions')
    #     print('Choose function: {}'.format(fun.__name__))
    #     fun(**d2).run()


if __name__ == '__main__':
    main()

# EOF