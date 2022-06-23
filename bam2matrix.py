#!/usr/bin/env python
#-*- encoding:utf-8 -*-
"""
Convert BAM to matrix
BAM --bam2bw()-- BigWig --bw2matrix()-- Matrix
"""


import os
import argparse
from bam2bw import Bam2bw
from bw2matrix import Bw2matrix
from utils import (make_config, update_obj, file_abspath, fix_out_dir, log,)
from parse_args import add_io_parser, add_bam_parser, add_bw_parser, get_init_args


class Bam2matrix(object):
    """
    Convert BAM to matrix:
    1. strand-specific: yes (prefix_sens.mat.gz, prefix_anti.mat.gz)
    2. strand-specific: no (prefix.mat.gz)
    """
    def __init__(self, **kwargs):
        # c = make_config(**kwargs)
        c = get_init_args(**kwargs)
        self = update_obj(self, c, force=True)
        self.init_args()


    def init_args(self):
        # bam_list, region_list, strand_specific
        self.out_dir = fix_out_dir(self.out_dir)
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.prefix)
        self.bam_list = file_abspath(self.bam_list)
        self.region_list = file_abspath(self.region_list)


    def run(self):
        # 1. bam to bigWig
        args = self.__dict__.copy()
        args.update({'out_dir': os.path.join(self.out_dir, '1.bam2bw')})
        bw_list = Bam2bw(**args).run()
        # 2. bigWig to matrix
        if self.strand_specific:
            args.update({
                'bw_list': None, 
                'bw_fwd_list': [i[0] for i in bw_list],
                'bw_rev_list': [i[1] for i in bw_list],
            })
        else:
            args.update({
                'bw_list': bw_list,
                'bw_fwd_list': None,
                'bw_rev_list': None,
            })
        args.update({'out_dir': os.path.join(self.out_dir, '2.bw2matrix')})
        return Bw2matrix(**args).run() # matrix, str or list


def get_args():
    example = ' '.join([
        'Example: \n',
        '$ python bam2matrix.py',
        '-b f1.bam f2.bam -r g1.bed -o out_dir -op metaplot',
        '-t scale-regions -u 2000 -d 2000 -m 2000 --binSize 100',
        '--blackListFileName bl.bed',
        '-sl f1 f2 -rl g1',
        '--startLabel TSS --endLabel TES -p 8',
    ])
    parser = argparse.ArgumentParser(
        prog='bam2matrix.py', description='bam2matrix', epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser = add_io_parser(parser)
    parser = add_bam_parser(parser)
    parser = add_bw_parser(parser)
    # parser.add_argument('-b', dest='bam_list', nargs='+', required=True,
    #     help='bam files')
    # parser.add_argument('-r', dest='region_list', nargs='+', required=True,
    #     help='region files, BED format')
    # parser.add_argument('-o', dest='out_dir', required=True,
    #     help='directory to output file')
    # parser.add_argument('-op', '--out-prefix', dest='prefix', default='bam2matrix',
    #     help='prefix for output files, default: [bam2matrix]')
    # parser.add_argument('-ss','--strand-specific', dest='strand_specific',
    #     action='store_true', help='Strand-specific, dUTP library')
    # parser.add_argument('-t', '--matrix-type', dest='matrix_type',
    #     default='scale-regions', choices=['scale-regions', 'reference-point'],
    #     help='choose the matrix type, default: [scale-regions]')
    # parser.add_argument('-bs', '--binSize', dest='binSize', type=int, default=50,
    #     help='the bin_size, default [50]')
    # parser.add_argument('-n', '--normalizeUsing', default='None', 
    #     choices=['RPKM', 'CPM', 'BPM', 'RPGC', 'None'],
    #     help='Use one of the method to normalize reads, default: ["None"]')
    # parser.add_argument('-g', '--genome', default=None,
    #     help='The reference genome of bam files, default [None]')
    # parser.add_argument('-es', '--effsize', dest='effectiveGenomeSize', type=int,
    #     default=None,
    #     help='effective genome size, if not specified, parse from bam header')
    # parser.add_argument('-sl', '--samplesLabel', nargs='+', default=None,
    #     help='labels for samples in plot, defautl: [None] auto')
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
    Bam2matrix(**args).run()


if __name__ == '__main__':
    main()


# EOF