#!/usr/bin/env python

"""
Convert BAM to profile
attention: strand-specific option
"""


import os
import sys
import pathlib
import argparse
import shutil

from utils import (
    make_config, update_obj, file_abspath, file_prefix, symlink_file,
    fix_label, fix_bw, is_valid_bam, is_valid_bigwig, is_valid_bed, is_valid_file,
    log, load_matrix_header
)
from bam2bw import Bam2bw
from bam2matrix import Bam2matrix
from bw2matrix import Bw2matrix
from matrix2profile import Matrix2profile


class Bam2profile(object):
    """
    Convert BAM to profile: strand-specific
    """
    def __init__(self, **kwargs):
        c = make_config(**kwargs)
        self = update_obj(self, c, force=True)
        self.init_args()
        self.init_files()
        self.check_files()
        
    
    def init_args(self):
        # strand specific
        if not hasattr(self, 'strand_specific'):
            self.strand_specific = False
        # check sample files (list)
        if is_valid_file(self.bam_list, is_valid_bam):
            if isinstance(self.bam_list, str):
                self.bam_list = [self.bam_list]
            self.bam_list = list(map(file_abspath, self.bam_list))
        else:
            raise ValueError('illegal, bam_list={}'.format(self.bam_list))
        # check region files (list)        
        if is_valid_file(self.region_list, is_valid_bed):
            if isinstance(self.region_list, str):
                self.region_list = [self.region_list]
            self.region_list = list(map(file_abspath, self.region_list))
        else:
            raise ValueError('illegal, region_list={}'.format(
                self.region_list
            ))
        # prefix
        if not isinstance(self.out_prefix, str):
            self.out_prefix = 'metplot'
            log.warning('out_prefix illegal, use "metaplot"')
        # perGroup
        self.perGroup = True # force
        # labels        
        self.samplesLabel = fix_label(self.bam_list, self.samplesLabel)
        self.regionsLabel = fix_label(self.region_list, self.regionsLabel)
    
    
    def check_files(self):
        # n_samples
        samples = self.bam_list
        n_samples = len(samples)
        if not len(samples) == len(self.samplesLabel):
            raise ValueError('samples: {} and lables: {} not identical'.format(
                len(samples), len(self.samplesLabel)
            ))
        msg = '\n'.join([
            '-'*80,
            '>> Sample files:',
            '\n'.join([
                '{} : {}'.format(i, k) for i,k in zip(self.samplesLabel, samples)
            ]),
            '-'*30,
            '>> Region files:',
            '\n'.join([
                '{} : {}'.format(i, k) for i,k in zip(self.regionsLabel, self.region_list)
            ]),
            '-'*80,
        ])
        return msg

    
    def init_files(self):
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.out_prefix)
        prefix = os.path.join(self.project_dir, self.out_prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.yaml'),
            'bed_dir': os.path.join(self.project_dir, '1.bed_files'),
            'bam_dir': os.path.join(self.project_dir, '2.bam_files'),
            'bw_dir':  os.path.join(self.project_dir, '3.bw_files'),
            'matrix_dir': os.path.join(self.project_dir, '4.matrix_files'),
            'profile_dir': os.path.join(self.project_dir, '5.profile_files'),
            'matrix_cmd': os.path.join(self.project_dir, prefix+'.matrix.sh'),
            'matrix': prefix+'.mat.gz',
            'matrix_sens': prefix+'.mat.sens.gz',
            'matrix_anti': prefix+'.mat.anti.gz',
            'profile_file': prefix+'.profile.pdf',
            'profile_file_sens': prefix+'.profile_sens.pdf',
            'profile_file_anti': prefix+'.profile_anti.pdf',
            'stdout': prefix+'.stdout',
            'stderr': prefix+'.stderr',
        }
        self = update_obj(self, args, force=True)
        for d in [self.bed_dir, self.bam_dir, self.bw_dir, self.matrix_dir, self.profile_dir]:
            if not os.path.exists(d):
                os.makedirs(d)   
        

    def bam2profile_ns(self):
        """
        strand-specific: False
        """
        # 1. bam to bw
        args1 = self.__dict__.copy()
        args1.update({
            'bam': self.bam_list,
            'out_dir': self.bw_dir,
            'out_prefix': None, # update for bw files
        })
        r1 = Bam2bw(**args1)
        bw_list = r1.run()
        
        # 2. bw to matrix
        args2 = self.__dict__.copy()
        args2.update({
            'bw_list': bw_list,
            'out_dir': self.matrix_dir,
        })
        r2 = Bw2matrix(**args2)
        r2.run()
        
        # 3. matrix to profile
        args3 = self.__dict__.copy()
        args3.update({
            'matrix': r2.matrix,
            'out_dir': self.profile_dir,
            'out_prefix': self.out_prefix,
        })
        r3 = Matrix2profile(**args3)
        r3.run()

        # 4. copy files
        symlink_file(r3.profile_file, self.profile_file)        
        return self.profile_file
        
        
    def bam2profile_ss(self):
        """
        strand-specific: True
        """
        # 1. bam to bw
        args1 = self.__dict__.copy()
        args1.update({
            'bam': self.bam_list,
            'out_dir': self.bw_dir,
            'out_prefix': None, # update for bw files
        })
        r1 = Bam2bw(**args1)
        bw_list = r1.run()
        bw_fwd_list = [i[0] for i in bw_list]
        bw_rev_list = [i[1] for i in bw_list]

        # 2. bw to matrix
        args2 = self.__dict__.copy()
        args2.update({
            'bw_list': None,
            'bw_fwd_list': bw_fwd_list,
            'bw_rev_list': bw_rev_list,
            'out_dir': self.matrix_dir,
        })
        r2 = Bw2matrix(**args2)
        r2.run()
        
        # 3. matrix to profile
        ## sens
        args3 = self.__dict__.copy()
        args3.update({
            'matrix': r2.matrix_sens,
            'out_dir': self.profile_dir,
            'out_prefix': self.out_prefix+'_sens',
        })
        r3 = Matrix2profile(**args3)
        r3.run()
        ## anti
        args4 = self.__dict__.copy()
        args4.update({
            'matrix': r2.matrix_anti,
            'out_dir': self.profile_dir,
            'out_prefix': self.out_prefix+'_anti',
        })
        r4 = Matrix2profile(**args4)
        r4.run()
        
        # 4. copy files
        symlink_file(r3.profile_file_sens, self.profile_file_sens)
        symlink_file(r3.profile_file_anti, self.profile_file_anti)        
        return [self.profile_file_sens, self.profile_file_anti]
        

    def run(self):
        fun = self.bam2profile_ss if self.strand_specific else self.bam2profile_ns
        return fun()


def get_args():
    example = ' '.join([
        'Example: \n',
        '$ python bam2profile.py', 
        '-b f1.bam f2.bam -r g1.bed -o out_dir --out-prefix metaplot', 
        '--matrix-type scale-regions -u 2000 -d 2000 -m 2000 --binSize 100',
        '--blackListFileName bl.bed', 
        '--samplesLabel f1 f2 --regionsLabel g1 g2',
        '--startLabel TSS --endLabel TES -p 8',
    ])
    parser = argparse.ArgumentParser(prog='bw2matrix.py',
                                     description='bw2matrix',
                                     epilog=example,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-b', dest='bam_list', nargs='+', required=True,
        help='bam files')
    parser.add_argument('-r', dest='region_list', nargs='+', required=True,
        help='region files')
    parser.add_argument('-o', dest='out_dir', required=True,
        help='directory to save files')
    parser.add_argument('--out-prefix', dest='out_prefix', default='bw2matrix',
        help='prefix for output files, default: [bw2matrix]')
    
    parser.add_argument('-t', '--matrix-type', dest='matrix_type', 
        default='scale-regions', choices=['scale-regions', 'reference-point'],
        help='choose the matrix type, default: [scale-regions]')
    parser.add_argument('-ss','--strand-specific', dest='strand_specific', 
        action='store_true', help='Strand-specific, dUTP library')
    parser.add_argument('-s', '--scaleFactor', dest='scaleFactor', type=float, 
        default=1.0, help='scale factor for the bam, default: [1.0]') 
    parser.add_argument('-n', '--normalizeUsing', default='None', 
        choices=['RPKM', 'CPM', 'BPM', 'RPGC', 'None'],
        help='Use one of the method to normalize reads, default: ["None"]')
    parser.add_argument('-es', '--effsize', dest='effectiveGenomeSize', type=int,
        default=None,
        help='effective genome size, if not specified, parse from bam header')
    parser.add_argument('-g', '--genome', default=None,
        help='The reference genome of bam files, default [None]')
    
    parser.add_argument('--samplesLabel', nargs='+', default=None,
        help='labels for samples in plot, defautl: [None] auto')
    parser.add_argument('--regionsLabel', nargs='+', default=None,
        help='labels for regions in plot, defautl: [None] auto')
    parser.add_argument('-bs', '--binSize', dest='binSize', type=int, default=50,
        help='the bin_size, default [50]')
    parser.add_argument('-u', '--beforeRegionStartLength', type=int, default=500,
        help='Distance upstream of TSS, default: [500]')
    parser.add_argument('-d', '--afterRegionStartLength', type=int, default=500,
        help='Distance downstream of TES, default: [500]')
    parser.add_argument('--startLabel', default='TSS',
        help='start label, default: [TSS]')
    parser.add_argument('--endLabel', default='TES',
        help='end label, default: [TES]')    
    parser.add_argument('-m', '--regionBodyLength', type=int, default=1000,
        help='Distance for all regions, default: [1000]')
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
    Bam2profile(**args).run()

    
if __name__ == '__main__':
    main()

# EOF