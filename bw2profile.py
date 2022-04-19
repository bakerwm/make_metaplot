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
    make_config, update_obj, load_yaml, dump_yaml, file_abspath, file_prefix,
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
        else:
            raise ValueError('illegal, region_list={}'.format(
                self.region_list))
        # bw files
        if is_valid_file(self.bw_fwd_list, is_valid_bigwig) and \
            is_valid_file(self.bw_rev_list, is_valid_bigwig):
            self.strand_specific = True
            tag = 0
            for f,r in zip(self.bw_fwd_list, self.bw_rev_list):
                f1 = file_prefix(f).strip('_fwd')
                r1 = file_prefix(r).strip('_rev')
                if not f1 == r1:
                    tag += 1
                    print('bw fwd/rev not matched: {} {}'.format(f, r))
            if tag > 0:
                raise ValueError('bw fwd/rev files illegal')
            # update samplesLabel
            self.samplesLabel = fix_label(self.bw_fwd_list, self.samplesLabel)
            self.samplesLabel = [i.rstrip('_fwd') for i in self.samplesLabel]
        elif is_valid_file(self.bw_list, is_valid_bigwig):
            self.strand_specific = False
            self.samplesLabel = fix_label(self.bw_rwd_list, self.samplesLabel)
        else:
            raise ValueError('bw_list, bw_fwd_list, bw_rev_list required:')
        # prefix
        if not isinstance(self.out_prefix, str):
            raise ValueError('out_prefix, expect str, got {}'.format(
                type(self.out_prefix).__name__))
        # labels        
        # self.samplesLabel = fix_label(self.bam_list, self.samplesLabel)
        self.regionsLabel = fix_label(self.region_list, self.regionsLabel)

            
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

            
    def split_bed(self, x, strand='+'):
        """
        Split BED file by strand
        * : original 
        + : _fwd
        - : _rev
        """
        xname = file_prefix(x)
        t = '_rev' if strand == '-' else '_fwd' if strand == '+' else ''
        out = os.path.join(self.bed_dir, xname+t+'.bed')
        if os.path.exists(out):
            print('file exists, skipped {}'.format(out))
        else:
            with open(x) as r, open(out, 'wt') as w:
                for line in r:
                    s = line.strip().split('\t')
                    if len(s) > 5:
                        if s[5] == strand or strand == '*':
                            w.write(line)
                    elif strand == '*':
                        w.write(line)
                    else:
                        pass
        return out

        
    def subset_matrix(self, x, o, samples=None, groups=None):
        """
        Example:
        computeMatrixOperations subset -m input.mat.gz -o output.mat.gz --groups "group 1"  --samples "sample 3"
        """
        d = load_matrix_header(x)
        args = ''
        if isinstance(samples, list):
            args += ' --samples {}'.format(' '.join(samples))
        if isinstance(groups, list):
            args += ' --groups {}'.format(' '.join(groups))
        if args == '': # empty
            if os.path.exists(o) and not self.overwrite:
                # if re-cal required, remove the old file
                log.info('subset_matrix() skipped, file exists: {}'.format(o))
            else:
                shutil.copy(x, o) # copy
        else:
            cmd = ' '.join([
                '{} subset'.format(shutil.which('computeMatrixOperations')),
                '-m {}'.format(x),
                '-o {}'.format(o),
                args
            ])
            # save command
            cmd_txt = os.path.join(os.path.dirname(o), file_prefix(o)+'.subset_matrix.sh')
            with open(cmd_txt, 'wt') as w:
                w.write(cmd+'\n')
            # run
            if os.path.exists(o) and not self.overwrite:
                # if re-cal required, remove the old file
                log.info('subset_matrix() skipped, file exists: {}'.format(o))
            else:
                log.info('run subset_matrix: {}'.format(o))
                os.system(cmd) # try-except ?

  
    def rbind_matrix(self, x, o):
        """
        Example:
        computeMatrixOperations rbind -m input1.mat.gz input2.mat.gz -o output.mat.gz
        """
        cmd = ' '.join([
            '{} rbind'.format(shutil.which('computeMatrixOperations')),
            '-m {}'.format(' '.join(x)),
            '-o {}'.format(o),
        ])
        # save command
        cmd_txt = os.path.join(os.path.dirname(o), file_prefix(o)+'.rbind_matrix.sh')
        with open(cmd_txt, 'wt') as w:
            w.write(cmd+'\n')
        # run
        if os.path.exists(o) and not self.overwrite:
            log.info('rbind_matrix() skipped, file exists: {}'.format(o))
        else:
            log.info('run rbind_matrix: {}'.format(o))
            os.system(cmd)

    
    def bw2matrix_ss(self):
        """
        Strand-specific: True
        """
        # 1. split regions file, by strand
        args1 = self.__dict__.copy()
        region_fwd_list = [self.split_bed(i, strand='+') for i in self.region_list]
        region_rev_list = [self.split_bed(i, strand='-') for i in self.region_list]
        
        # labels
        sl_fwd = [i+'_fwd' for i in self.samplesLabel]
        sl_rev = [i+'_rev' for i in self.samplesLabel]
        rl_fwd = [os.path.basename(i) for i in region_fwd_list] # 
        rl_rev = [os.path.basename(i) for i in region_rev_list] #

        # 2. bw to matrix
        args2 = self.__dict__.copy()
        args2.update({
            'out_dir': self.matrix_dir,
            'out_prefix': self.out_prefix,
            'bw_list': self.bw_fwd_list + self.bw_rev_list,
            'region_list': region_fwd_list + region_rev_list,
            'samplesLabel': sl_fwd + sl_rev,
            'regionsLabel': rl_fwd + rl_rev,
        })
        r2 = Bw2matrix(**args2)
        r2.run()
        
        # 3. subset matrix
        ## sense
        sens1 = os.path.join(args2['out_dir'], self.out_prefix+'.sens_1.mat.gz')
        sens2 = os.path.join(args2['out_dir'], self.out_prefix+'.sens_2.mat.gz')
        self.subset_matrix(r2.matrix, sens1, samples=sl_fwd, groups=rl_fwd) # fwd+fwd
        self.subset_matrix(r2.matrix, sens2, samples=sl_rev, groups=rl_rev) # rev+rev
        ## antisense
        anti1 = os.path.join(args2['out_dir'], self.out_prefix+'.anti_1.mat.gz')
        anti2 = os.path.join(args2['out_dir'], self.out_prefix+'.anti_2.mat.gz')
        self.subset_matrix(r2.matrix, anti1, samples=sl_fwd, groups=rl_rev) # fwd+rev
        self.subset_matrix(r2.matrix, anti2, samples=sl_rev, groups=rl_fwd) # rev+fwd
        # 4. merge
#         self.rbind_matrix(x=[sens1, sens2], o=self.matrix_sens)
#         self.rbind_matrix(x=[anti1, anti2], o=self.matrix_anti)
        Matrix_rbind(m=[sens1, sens2], o=self.matrix_sens).run()
        Matrix_rbind(m=[anti1, anti2], o=self.matrix_anti).run()
        return (self.matrix_sens, self.matrix_anti)
    

    def bw2matrix_ns(self):
        """
        Strand-specific: False
        """
        # 1. split regions file, by strand
        ## skipped
        
        # 2. bw to matrix
        args2 = self.__dict__.copy()
        args2.update({
            'out_dir': self.matrix_dir,
            'out_prefix': self.out_prefix,
            'bw_list': self.bw_list,
            'region_list': self.region_list,
            'samplesLabel': self.samplesLabel,
            'regionsLabel': self.regionsLabel,
        })
        r2 = Bw2matrix(**args2)
        r2.run()
        
        # 3. copy matrix file
        src_dir_rel = os.path.relpath(os.path.dirname(r2.matrix), os.path.dirname(self.matrix))        
        src_rel = os.path.join(src_dir_rel, os.path.basename(r2.matrix))
        if not os.path.exists(self.matrix):
            os.symlink(src_rel, self.matrix)
        
        return r2.matrix

    
    def run(self):
        if self.strand_specific:
            # 1. bw to matrix
            matrix_sens, matrix_anti = self.bw2matrix_ss()
            
            # 2. matrix to profile: sens
            args2 = self.__dict__.copy()
            args2.update({
                'matrix': matrix_sens,
                'out_dir': self.profile_dir,
                'perGroup': True,
            })
            r2 = Matrix2profile(**args2)
            r2.run()
            
            # 3. copy files: sens
            if not os.path.exists(self.profile_file_sens) or self.overwrite:
                shutil.copy(r2.profile_file, self.profile_file_sens)

            # 4 matrix to profile: anti
            args3 = self.__dict__.copy()
            args3.update({
                'matrix': matrix_anti,
                'out_dir': self.profile_dir,
                'perGroup': True,
            })
            r3 = Matrix2profile(**args3)
            r3.run()
            
            # 5 copy files: anti
            if not os.path.exists(self.profile_file_anti) or self.overwrite:
                shutil.copy(r3.profile_file, self.profile_file_anti)
        else:
            # 1. bw to matrix
            matrix = self.bw2matrix_ns()
            
            # 2. matrix to profile
            args2 = self.__dict__.copy()
            args2.update({
                'matrix': matrix,
                'out_dir': self.profile_dir,
                'perGroup': True,
            })
            r2 = Matrix2profile(**args2)
            r2.run()
            
            # 3. copy files
            if not os.path.exists(self.profile_file) or self.overwrite:
                shutil.copy(r2.profile_file, self.profile_file)


def get_args():
    example = ' '.join([
        '$ python bw2profile.py', 
        '-b f1.bw f2.bw -r g1.bed g2.bed -o out_dir --out-prefix metaplot', 
        '--matrix-type scale-regions -u 2000 -d 2000 -m 2000 --binSize 100',
        '--blackListFileName bl.bed', 
        '--samplesLabel f1 f2 --regionsLabel g1 g2',
        '--startLabel TSS --endLabel TES -p 8',
    ])
    parser = argparse.ArgumentParser(prog='bw2matrix.py',
                                     description='bw2matrix',
                                     epilog=example,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-b', dest='bw_list', nargs='+', required=False,
        help='bw files')
    parser.add_argument('-bf', dest='bw_fwd_list', nargs='+',
        help='bw files on forward strand')
    parser.add_argument('-br', dest='bw_rev_list', nargs='+',
        help='bw files on reverse strand')
    parser.add_argument('-r', dest='region_list', nargs='+', required=True,
        help='region files')
    parser.add_argument('-o', dest='out_dir', required=True,
        help='directory to save files')
    parser.add_argument('--out-prefix', dest='out_prefix', default='bw2matrix',
        help='prefix for output files, default: [bw2matrix]')
    parser.add_argument('--strand-specific', dest='strand_specific', 
        action='store_true', help='for strand-specific anslysis')
    parser.add_argument('-t', '--matrix-type', dest='matrix_type', 
        default='scale-regions', choices=['scale-regions', 'reference-point'],
        help='choose the matrix type, default: [scale-regions]')     
    parser.add_argument('--samplesLabel', nargs='+', default=None,
        help='labels for samples in plot, defautl: [None] auto')
    parser.add_argument('--regionsLabel', nargs='+', default=None,
        help='labels for regions in plot, defautl: [None] auto')
    parser.add_argument('--binSize', type=int, default=50,
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
    Bw2profile(**args).run()

    
if __name__ == '__main__':
    main()

# EOF