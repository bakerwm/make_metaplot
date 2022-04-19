#!/usr/bin/env python

"""
Convert BAM to matrix
attention: strand-specific option

bam2bw + bw2matrix
"""


import os
# import sys
import pathlib
import argparse
import shutil
import pysam

from utils import (
    make_config, update_obj, file_abspath, file_prefix, symlink_file,
    fix_label, fix_bw, is_valid_bam, is_valid_bigwig, is_valid_bed, is_valid_file,
    log, load_matrix_header
)
from bam2bw import Bam2bw
from bw2matrix import Bw2matrix
from matrix_rbind import Matrix_rbind


############################################################################################
# wrap all in one
#
# 1. bam -> plots (DNA + RNA)
# 2. bw -> plots (DNA)
class Bam2matrix(object):
    """
    Generate matrix for BAM file and regions file

    Arguments
    ---------
    bam_list:

    region_list:

    out_dir:

    out_prefix:

    Optional
    --------
    samplesLabel:
    regionsLabel:
    ...

    1. non-stranded: Bam2bw + Bw2matrix
    2. stranded: Bam2bw (sens + anti) + Bw2matrix (sens + anti) + mergeMatrix (rbind)
    """
    def __init__(self, **kwargs):
        c = make_config(**kwargs)
        self = update_obj(self, c, force=True)
        self.init_args()
        self.init_files()


    def init_args(self):
        if not hasattr(self, 'strand_specific'):
            self.strand_specific = False
        # region files (list)
        if is_valid_file(self.region_list, is_valid_bed):
            if isinstance(self.region_list, str):
                self.region_list = [self.region_list]
            self.region_list = list(map(file_abspath, self.region_list))
            self.regionsLabel = fix_label(self.region_list, self.regionsLabel)
        else:
            raise ValueError('illegal, region_list={}'.format(
                self.region_list
            ))
        # score files (list)
        if is_valid_file(self.bam_list, is_valid_bam):
            if isinstance(self.bam_list, str):
                self.bam_list = [self.bam_list]
            self.bam_list = list(map(file_abspath, self.bam_list))
        else:
            raise ValueError('illegal, bam_list={}'.format(self.bam_list))
        # prefix
        if not isinstance(self.out_prefix, str):
            self.out_prefix = 'metaplot'
            log.warning('out_prefix failed, use "metaplot" instead')
        # labels
        self.samplesLabel = fix_label(self.bam_list, self.samplesLabel)
        self.regionsLabel = fix_label(self.region_list, self.regionsLabel)
        # type
        if not self.matrix_type in ['scale-regions', 'reference-point']:
            self.matrix_type = 'scale-regions'


    def init_files(self):
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.out_prefix)
        prefix = os.path.join(self.project_dir, self.out_prefix)
        args = {
            'bw_dir': os.path.join(self.project_dir, '1.bw_files'),
            'matrix_dir': os.path.join(self.project_dir, '2.matrix_files'),
            'matrix': prefix+'.mat.gz',
            'matrix_sens': prefix+'_sens.mat.gz',
            'matrix_anti': prefix+'_anti.mat.gz',
            'stdout': prefix+'.matrix.stdout',
            'stderr': prefix+'.matrix.stderr',
        }
        self = update_obj(self, args, force=True)
        for d in [self.bw_dir, self.matrix_dir]:
            if not os.path.exists(d):
                os.makedirs(d)


    def check_files(self):
        msg = '\n'.join([
            '-'*80,
            '>> BAM files:',
            '\n'.join([
                '{} : {}'.format(i, k) for i,k in zip(self.samplesLabel, self.bam_list)
            ]),
            '-'*30,
            '>> region files:',
            '\n'.join([
                '{} : {}'.format(i, k) for i,k in zip(self.regionsLabel, self.region_list)
            ]),
            '-'*80,
        ])
        return msg


    def bam2matrix_ns(self):
        """
        strand-specific: False
        """
        # 1. bam2bw
        args1 = self.__dict__.copy()
        args1.update({
            'bam': self.bam_list,
            'out_dir': self.bw_dir,
        })
        r1 = Bam2bw(**args1)
        bw_list = r1.run() # bw_list

        # 2. bw2matrix
        args2 = self.__dict__.copy()
        args2.update({
            'bw_list': bw_list,
            'out_dir': self.matrix_dir,
        })
        r2 = Bw2matrix(**args2)
        matrix = r2.run()

        # 3. copy files
        symlink_file(matrix, self.matrix)

        # 4. check file
        if not os.path.exists(self.matrix):
            log.error('bam2matrix() failed, file not exists: {}'.format(
                self.matrix
            ))
        return matrix # str


    def bam2matrix_ss(self):
        """
        strand-specific: True
        """
        # 1. bam2bw
        args1 = self.__dict__.copy()
        args1.update({
            'bam': self.bam_list,
            'out_dir': self.bw_dir,
            'out_prefix': None, # auto
        })
        r1 = Bam2bw(**args1)
        bw_list = r1.run() # bw_list

        # 2. bw2matrix
        args2 = self.__dict__.copy()
        args2.update({
            'bw_list': None,
            'bw_fwd_list': [i[0] for i in bw_list],
            'bw_rev_list': [i[1] for i in bw_list],
            'out_dir': self.matrix_dir,
        })
        r2 = Bw2matrix(**args2)
        matrix = r2.run() # tuple? why
        matrix = list(matrix) 


        # 3. copy files
        symlink_file(matrix[0], self.matrix_sens)
        symlink_file(matrix[1], self.matrix_anti)

        # 4. check file
        if not is_valid_file(matrix, os.path.exists):
            log.error('bam2matrix() failed, file not exists: {}, {}'.format(
                self.matrix_sens, self.matrix_anti
            ))
        return matrix # list, [sens, anti]


    def run(self):
        print(self.check_files())
        fun = self.bam2matrix_ss if self.strand_specific else self.bam2matrix_ns
        return fun()


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
    parser.add_argument('-b', dest='bam_list', nargs='+', required=True,
        help='bam files')
    parser.add_argument('-r', dest='region_list', nargs='+', required=True,
        help='region files, BED format')
    parser.add_argument('-o', dest='out_dir', required=True,
        help='directory to output file')
    parser.add_argument('-op', '--out-prefix', dest='out_prefix', default='bam2matrix',
        help='prefix for output files, default: [bam2matrix]')
    parser.add_argument('-ss','--strand-specific', dest='strand_specific',
        action='store_true', help='Strand-specific, dUTP library')
    parser.add_argument('-t', '--matrix-type', dest='matrix_type',
        default='scale-regions', choices=['scale-regions', 'reference-point'],
        help='choose the matrix type, default: [scale-regions]')
    parser.add_argument('-bs', '--binSize', dest='binSize', type=int, default=50,
        help='the bin_size, default [50]')
    parser.add_argument('-n', '--normalizeUsing', default='None', 
        choices=['RPKM', 'CPM', 'BPM', 'RPGC', 'None'],
        help='Use one of the method to normalize reads, default: ["None"]')
    parser.add_argument('-g', '--genome', default=None,
        help='The reference genome of bam files, default [None]')
    parser.add_argument('-es', '--effsize', dest='effectiveGenomeSize', type=int,
        default=None,
        help='effective genome size, if not specified, parse from bam header')
    parser.add_argument('-sl', '--samplesLabel', nargs='+', default=None,
        help='labels for samples in plot, defautl: [None] auto')
    parser.add_argument('-u', '--beforeRegionStartLength', type=int, default=500,
        help='Distance upstream of TSS, default: [500]')
    parser.add_argument('-d', '--afterRegionStartLength', type=int, default=500,
        help='Distance downstream of TES, default: [500]')
    parser.add_argument('-m', '--regionBodyLength', type=int, default=1000,
        help='Distance for all regions, default: [1000]')
    parser.add_argument('-st', '--startLabel', default='TSS',
        help='start label, default: [TSS]')
    parser.add_argument('-ed', '--endLabel', default='TES',
        help='end label, default: [TES]')
    parser.add_argument('--sortRegions', default='keep',
        choices=['descend', 'ascend', 'no', 'keep'],
        help='The output should be sorted by the way.')
    parser.add_argument('--sortUsing', default='mean',
        choices=['mean', 'median', 'max', 'min', 'sum', 'region_length'],
        help='Which method should be used for sorting, default: [mean]')
    parser.add_argument('--sortUsingSamples', type=str, default=None,
        help='List of sample numbers for sorting, default: [None]')
    parser.add_argument('--averageTypeBins',
        choices=['mean', 'median', 'min', 'max', 'std', 'sum'],
        help='Define the type of method for bin size range, default: [mean]')
    parser.add_argument('-bl', '--blackListFileName', dest='blackListFileName',
        default=None, help='blacklist file')
    parser.add_argument('-p', dest='numberOfProcessors', type=int, default=4,
        help='number of processors, default: [4]')
    parser.add_argument('-O', '--overwrite', dest='overlap', action='store_true',
        help='Overwrite output file')
    return parser


def main():
    args = vars(get_args().parse_args())
    Bam2matrix(**args).run()


if __name__ == '__main__':
    main()

# EOF