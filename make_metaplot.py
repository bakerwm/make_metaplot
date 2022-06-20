#!/usr/bin/env python

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
   matrix-1 + matrix-2 -> matrix_fwd
   matrix-3 + matrix-4 -> matrix_rev
7. matrix (fwd/rev) -> plots

## Examples

### A. Simple example for DNA
(bam/bw files; region files; out_dir)


### B. Simple example for RNA
(bam files(required); region files; out_dir)


### C. Complicated examples
(normalization, scale, ...)
"""

import os
# import sys
# import pathlib
import argparse
# import shutil
# import logging
# from matplotlib import colors
# from xopen import xopen
# import pysam
# import pyBigWig
# import json
# import yaml
# import toml

from bam2bw import Bam2bw
from bam2matrix import Bam2matrix
from bw2matrix import Bw2matrix
from bam2profile import Bam2profile
from bw2profile import Bw2profile
from matrix2profile import Matrix2profile
from matrix2heatmap import Matrix2heatmap
from utils import (
    make_config, update_obj, load_yaml, dump_yaml, file_abspath, file_prefix,
    fix_label, fix_bw, is_valid_bam, is_valid_bigwig, is_valid_matrix, log,
    is_valid_file
)


"""
# Examples
$ computeMatrixOperations subset -m foo.mat.gz -o forward.mat.gz --samples SRR648667.forward SRR648668.forward SRR648669.forward SRR648670.forward
$ computeMatrixOperations subset -m foo.mat.gz -o reverse.mat.gz --samples SRR648667.reverse SRR648668.reverse SRR648669.reverse SRR648670.reverse

# Examples
$ computeMatrixOperations filterStrand -m forward.mat.gz -o forward.subset.mat.gz --strand -
$ computeMatrixOperations filterStrand -m reverse.mat.gz -o reverse.subset.mat.gz --strand +

# Examples
$ computeMatrixOperations rbind -m forward.subset.mat.gz reverse.subset.mat.gz -o merged.mat.gz
$ computeMatrixOperations sort -m merged.mat.gz -o sorted.mat.gz -R genes.gtf
"""


def find_fun(x):
    """
    Generate metaplot file

    Parameters
    ----------
    x : dict
        The config, in dict

    from x to profile.
    Find the proper function to gene metaplot(profile) file

    priority:
    1. Matrix2profile: matrix
    2. Bw2profile: bw_fwd_list, bw_rev_list
    3. Bw2profile: bw_list
    4. Bam2profile: bam_list
    """
#     try:
#         args = load_yaml(x)
#     except:
#         log.error('unable to load YAML file: {}'.format(x))
#         return None
    args = x #
    # variables
    m = args.get('matrix', None)
    bf = args.get('bw_fwd_list', None)
    br = args.get('bw_rev_list', None)
    bw = args.get('bw_list', None)
    bam = args.get('bam_list', None)
    ## choose
    # 1. matrix
    if is_valid_matrix(m):
        fun = Matrix2profile
        args.update({
            'bw_fwd_list': None,
            'bw_rev_list': None,
            'bw_list': None,
            'bam_list': None,
        })
    elif is_valid_file(bf, is_valid_bigwig) and is_valid_file(br, is_valid_bigwig):
        fun = Bw2profile
        args.update({
            'matrix': None,
            'bw_list': None,
            'bam_list': None,
        })
    elif is_valid_file(bw, is_valid_bigwig):
        fun = Bw2profile
        args.update({
            'matrix': None,
            'bw_fwd_list': None,
            'bw_rev_list': None,
            'bam_list': None,
        })
    elif is_valid_file(bam, is_valid_bam):
        fun = Bam2profile
        args.update({
            'matrix': None,
            'bw_fwd_list': None,
            'bw_rev_list': None,
            'bw_list': None,
        })
    else:
        fun = None
    return (fun, args)


def show_help(x):
    msg = '\n'.join([
        '-'*80,
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
        '  - yMin, yMax',
        '# 3. Run the command again',
        '$ python make_profile.py -c {}'.format(x),
        '-'*80,
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
    args = {
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
            dump_yaml(make_config(), args['config'])
        show_help(args['config'])
    else:
        d1 = make_config()
        if d1 is None:
            print('error, check config file: {}'.format(args['config']))
            sys.exit(1)
        d1.update(load_yaml(args['config']))
        fun, d2 = find_fun(d1)
        if fun is None:
            raise ValueError('unable to determin functions')
        print('Choose function: {}'.format(fun.__name__))
        fun(**d2).run()


if __name__ == '__main__':
    main()

# EOF