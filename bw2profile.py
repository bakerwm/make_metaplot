#!/usr/bin/env python

"""
Convert BW to profile
attention: strand-specific option
- bw_fwd_list
- bw_rev_list
- bw_list
"""


import os
import pathlib
import argparse
import shutil

from utils import (
    make_config, update_obj, file_abspath, file_prefix, symlink_file,
    fix_label, fix_bw, is_valid_bam, is_valid_bigwig, is_valid_bed, is_valid_file,
    log, load_matrix_header
)
# from bam2bw import Bam2bw
from bw2matrix import Bw2matrix
from matrix2profile import Matrix2profile
from matrix_rbind import Matrix_rbind


class Bw2profile(object):
    """
    Convert BigWig to profile
    require: bw_fwd_list, bw_rev_list, bw_list
    """    
    def __init__(self, **kwargs):
        c = make_config(**kwargs)
        self = update_obj(self, c, force=True)
        self.init_args()
        self.init_files()


    def init_args(self):
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
        if is_valid_file(self.bw_fwd_list, is_valid_bigwig) and \
            is_valid_file(self.bw_rev_list, is_valid_bigwig):
            self.strand_specific = True
            tag = 0
            for f,r in zip(self.bw_fwd_list, self.bw_rev_list):
                f1 = file_prefix(f).replace('_fwd', '')
                r1 = file_prefix(r).replace('_rev', '')
                if not f1 == r1:
                    tag += 1
                    msg = '\n'.join([
                        '{} : {}'.format(f1, f),
                        '{} : {}'.format(r1, r),
                    ])
                    print(msg)
                    print('bw fwd/rev not matched')
            if tag > 0:
                msg = '\n'.join([
                    'fwd : {}'.format(';'.join(f)),
                    'rev : {}'.format(';'.join(r)),
                ])
                print(msg)
                raise ValueError('bw fwd/rev files illegal')
            # update samplesLabel
            self.samplesLabel = fix_label(self.bw_fwd_list, self.samplesLabel)
            self.samplesLabel = [i.rstrip('_fwd') for i in self.samplesLabel]
        elif is_valid_file(self.bw_list, is_valid_bigwig):
            self.strand_specific = False
            self.samplesLabel = fix_label(self.bw_list, self.samplesLabel)
        else:
            raise ValueError('bw_list, bw_fwd_list, bw_rev_list required:')
        # prefix
        if not isinstance(self.out_prefix, str):
            self.out_prefix = 'metaplot'
            log.warning('out_prefix failed, use "metaplot" instead')
        # labels
        self.regionsLabel = fix_label(self.region_list, self.regionsLabel)
        # type
        if not self.matrix_type in ['scale-regions', 'reference-point']:
            self.matrix_type = 'scale-regions'
        # perGroup
        self.perGroup = True # force

            
    def init_files(self):
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.out_prefix)
        prefix = os.path.join(self.project_dir, self.out_prefix)
        args = {
            # 'config': os.path.join(self.project_dir, 'config.yaml'),
            # 'bed_dir': os.path.join(self.project_dir, '1.bed_files'),
            # 'bam_dir': os.path.join(self.project_dir, '2.bam_files'),
            # 'bw_dir':  os.path.join(self.project_dir, '3.bw_files'),
            'matrix_dir': os.path.join(self.project_dir, '2.matrix_files'),
            'profile_dir': os.path.join(self.project_dir, '3.profile_files'),
            # 'matrix_cmd': os.path.join(self.project_dir, prefix+'.matrix.sh'),
            # 'matrix': prefix+'.mat.gz',
            # 'matrix_sens': prefix+'.mat.sens.gz',
            # 'matrix_anti': prefix+'.mat.anti.gz',
            'profile_file': prefix+'.profile.pdf',
            'profile_file_sens': prefix+'.profile_sens.pdf',
            'profile_file_anti': prefix+'.profile_anti.pdf',
            'stdout': prefix+'.stdout',
            'stderr': prefix+'.stderr',
        }
        self = update_obj(self, args, force=True)
        for d in [self.matrix_dir, self.profile_dir]:
            if not os.path.exists(d):
                os.makedirs(d)
                

    def check_files(self):
        msg = '\n'.join([
            '-'*80,
            '>> score files:',
            '\n'.join([
                '{} : {}'.format(i, k) for i,k in zip(self.samplesLabel, self.bw_list)
            ]),
            '-'*30,
            '>> region files:',
            '\n'.join([
                '{} : {}'.format(i, k) for i,k in zip(self.regionsLabel, self.region_list)
            ]),
            '-'*80,
        ])
        return msg
    
    
    def bw2profile_ns(self):
        """
        strand-specific: False
        """
        # 1. bw to matrix
        args1 = self.__dict__.copy()
        args1.update({
            'bw_list': self.bw_list,
            'out_dir': self.matrix_dir,
        })
        r1 = Bw2matrix(**args1)
        matrix = r1.run() # str

        # 2. matrix to profile
        args2 = self.__dict__.copy()
        args2.update({
            'matrix': matrix, # str
            'out_dir': self.profile_dir,
            # 'out_prefix': self.out_prefix,
        })
        r2 = Matrix2profile(**args2)
        r2.run()

        # 4. copy files
        symlink_file(r2.profile_file, self.profile_file)        
        return self.profile_file
        
        
    def bw2profile_ss(self):
        """
        strand-specific: True
        """
        # 1. bw to matrix
        args1 = self.__dict__.copy()
        args1.update({
            'bw_list': self.bw_list,
            'out_dir': self.matrix_dir,
        })
        r1 = Bw2matrix(**args1)
        matrix = r1.run() # tuple
        matrix = list(matrix)

        # 2. matrix to profile
        ## sens
        args2 = self.__dict__.copy()
        args2.update({
            'matrix': matrix[0],
            'out_dir': self.profile_dir,
            'out_prefix': self.out_prefix+'_sens',
        })
        r2 = Matrix2profile(**args2)
        r2.run()
        ## anti
        args3 = self.__dict__.copy()
        args3.update({
            'matrix': matrix[1],
            'out_dir': self.profile_dir,
            'out_prefix': self.out_prefix+'_anti',
        })
        r3 = Matrix2profile(**args3)
        r3.run()

        # 3. copy files
        symlink_file(r2.profile_file, self.profile_file_sens)
        symlink_file(r3.profile_file, self.profile_file_anti)        
        return [self.profile_file_sens, self.profile_file_anti]
    
    
    def run(self):
        fun = self.bw2profile_ss if self.strand_specific else self.bw2profile_ns
        return fun()


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
    parser.add_argument('-b', dest='bw_list', nargs='+', required=False,
        help='bw files')
    parser.add_argument('-bf', dest='bw_fwd_list', nargs='+', required=False,
        help='bw files on forward strand')
    parser.add_argument('-br', dest='bw_rev_list', nargs='+', required=False,
        help='bw files on reverse strand')
    parser.add_argument('-r', dest='region_list', nargs='+', required=True,
        help='region files')
    parser.add_argument('-o', dest='out_dir', required=True,
        help='directory to save bigWig file')
    parser.add_argument('-op', '--out-prefix', dest='out_prefix', default='bw2matrix',
        help='prefix for output files, default: [bw2matrix]')
    parser.add_argument('--plotType', default='lines',
        choices=['lines', 'fill', 'se', 'std', 'overlapped_lines', 'heatmap'],
        help='type of the plot, default: [lines]')
    parser.add_argument('--colors', nargs='+', default=None,
        help='colors for the lines, default: [None] auto')
    parser.add_argument('--averageType', default='mean',
        choices=['mean', 'median', 'max', 'min', 'sum', 'region_length'],
        help='Which method should be used for sorting, default: [mean]')    
    parser.add_argument('-t', '--matrix-type', dest='matrix_type',
        default='scale-regions', choices=['scales-regions', 'reference-point'],
        help='choose the matrix type, default: [scale-regions]')
    parser.add_argument('-sl', '--samplesLabel', nargs='+', default=None,
        help='samples label')
    parser.add_argument('-rl', '--regionsLabel', nargs='+', default=None,
        help='labels for regions in plot, defautl: [None] auto')
    parser.add_argument('-bs', '--binSize', type=int, default=50,
        help='the bin_size, default [50]')
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
    Bw2profile(**args).run()

    
if __name__ == '__main__':
    main()

# EOF