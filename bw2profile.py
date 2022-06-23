#!/usr/bin/env python
#-*- encoding:utf8 -*-
"""
Convert bigWig to profie
1. Bw2matrix() : bigWig to matrix
2. Matrix2profile() : matrix to profile
"""


import os
import argparse
from bw2matrix import Bw2matrix
from matrix2profile import Matrix2profile
from utils import (
    make_config, update_obj, fix_out_dir, log,
)
from parse_args import (
    add_io_parser, add_bw_parser, add_plot_parser, get_init_args,
)


class Bw2profile(object):
    """
    Convert BigWig to profile
    require: bw_fwd_list, bw_rev_list, bw_list, region_list
    """    
    def __init__(self, **kwargs):
        # c = make_config(**kwargs)
        c = get_init_args(**kwargs)
        self = update_obj(self, c, force=True)
        self.out_dir = fix_out_dir(self.out_dir)


    def run(self):
        # 1. bigWig to matrix
        args = self.__dict__.copy()
        args.update({
            'out_dir': os.path.join(self.out_dir, '2.bw2matrix'),
        })
        m = Bw2matrix(**args).run()
        # 2. matrix to profile
        args.update({
            'out_dir': os.path.join(self.out_dir, '3.matrix2profile'),
        })
        if isinstance(m, list):
            # sense strand
            args.update({
                'matrix': m[0],
                'prefix': None,
            })
            Matrix2profile(**args).run()
            # antisense strand
            args.update({'matrix': m[1]})
            Matrix2profile(**args).run()
        else:
            args.update({'matrix': m})
            Matrix2profile(**args).run()


def get_args():
    example = ' '.join([
        '$ python bw2profile.py', 
        '-b f1.bw f2.bw -r g1.bed g2.bed -o out_dir --out-prefix metaplot', 
        '--matrix-type scale-regions -u 2000 -d 2000 -m 2000 --binSize 100',
        '--blackListFileName bl.bed', 
        '--samplesLabel f1 f2 --regionsLabel g1 g2',
        '--startLabel TSS --endLabel TES -p 8',
    ])
    parser = argparse.ArgumentParser(
        prog='bw2matrix.py', description='bw2matrix', epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser = add_io_parser(parser)
    parser = add_bw_parser(parser)
    parser = add_plot_parser(parser)

    # parser.add_argument('-b', dest='bw_list', nargs='+', required=False,
    #     help='bw files')
    # parser.add_argument('-bf', dest='bw_fwd_list', nargs='+', required=False,
    #     help='bw files on forward strand')
    # parser.add_argument('-br', dest='bw_rev_list', nargs='+', required=False,
    #     help='bw files on reverse strand')
    # parser.add_argument('-r', dest='region_list', nargs='+', required=True,
    #     help='region files')
    # parser.add_argument('-o', dest='out_dir', required=True,
    #     help='directory to save bigWig file')
    # parser.add_argument('-op', '--out-prefix', dest='prefix', default='bw2matrix',
    #     help='prefix for output files, default: [bw2matrix]')
    # parser.add_argument('--plotType', default='lines',
    #     choices=['lines', 'fill', 'se', 'std', 'overlapped_lines', 'heatmap'],
    #     help='type of the plot, default: [lines]')
    # parser.add_argument('--colors', nargs='+', default=None,
    #     help='colors for the lines, default: [None] auto')
    # parser.add_argument('--averageType', default='mean',
    #     choices=['mean', 'median', 'max', 'min', 'sum', 'region_length'],
    #     help='Which method should be used for sorting, default: [mean]')    
    # parser.add_argument('-t', '--matrix-type', dest='matrix_type',
    #     default='scale-regions', choices=['scale-regions', 'reference-point'],
    #     help='choose the matrix type, default: [scale-regions]')
    # parser.add_argument('-sl', '--samplesLabel', nargs='+', default=None,
    #     help='samples label')
    # parser.add_argument('-rl', '--regionsLabel', nargs='+', default=None,
    #     help='labels for regions in plot, defautl: [None] auto')
    # parser.add_argument('-bs', '--binSize', type=int, default=50,
    #     help='the bin_size, default [50]')
    # parser.add_argument('-u', '--beforeRegionStartLength', type=int, default=500,
    #     help='Distance upstream of TSS, default: [500]')
    # parser.add_argument('-d', '--afterRegionStartLength', type=int, default=500,
    #     help='Distance downstream of TES, default: [500]')
    # parser.add_argument('-m', '--regionBodyLength', type=int, default=1000,
    #     help='Distance for all regions, default: [1000]')
    # parser.add_argument('-st', '--startLabel', default='TSS',
    #     help='start label, default: [TSS]')
    # parser.add_argument('-ed', '--endLabel', default='TES',
    #     help='end label, default: [TES]')
    # parser.add_argument('--sortRegions', default='keep',
    #     choices=['descend', 'ascend', 'no', 'keep'],
    #     help='The output should be sorted by the way.')
    # parser.add_argument('--sortUsing', default='mean',
    #     choices=['mean', 'median', 'max', 'min', 'sum', 'region_length'],
    #     help='Which method should be used for sorting, default: [mean]')
    # parser.add_argument('--sortUsingSamples', type=str, default=None,
    #     help='List of sample numbers for sorting, default: [None]')
    # parser.add_argument('--averageTypeBins',
    #     choices=['mean', 'median', 'min', 'max', 'std', 'sum'],
    #     help='Define the type of method for bin size range, default: [mean]')
    # parser.add_argument('-bl', '--blackListFileName', dest='blackListFileName',
    #     default=None, help='blacklist file')
    # parser.add_argument('-p', dest='numberOfProcessors', type=int, default=4,
    #     help='number of processors, default: [4]')
    # parser.add_argument('-O', '--overwrite', dest='overlap', action='store_true',
    #     help='Overwrite output file')
    return parser


def main():
    args = vars(get_args().parse_args())
    Bw2profile(**args).run()

    
if __name__ == '__main__':
    main()

# EOF