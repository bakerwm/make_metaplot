#!/usr/bin/env python

"""
Convert BAM to matrix
attention: strand-specific option
"""


import os
# import sys
import pathlib
import argparse
import shutil
import pysam

from utils import (
    make_config, update_obj, load_yaml, dump_yaml, file_abspath, file_prefix,
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
            if isinstance(self.bam_list, str):
                self.bam_list = [self.bam_list]
            self.region_list = list(map(file_abspath, self.region_list))
        else:
            raise ValueError('illegal, region_list={}'.format(
                self.region_list))
        # prefix
        if not isinstance(self.out_prefix, str):
            raise ValueError('out_prefix, expect str, got {}'.format(
                type(self.out_prefix).__name__))
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
            'config': os.path.join(self.project_dir, 'config.matrix.yaml'),
            'bed_dir': os.path.join(self.project_dir, '1.bed_files'),
            'bam_dir': os.path.join(self.project_dir, '2.bam_files'),
            'bw_dir':  os.path.join(self.project_dir, '3.bw_files'),
            'matrix_dir': os.path.join(self.project_dir, '4.matrix_files'),
            'matrix_cmd': os.path.join(self.project_dir, prefix+'.matrix.sh'),
            'matrix': prefix+'.mat.gz',
            'matrix_sens': prefix+'.mat.sens.gz',
            'matrix_anti': prefix+'.mat.anti.gz',
            'stdout': prefix+'.matrix.stdout',
            'stderr': prefix+'.matrix.stderr',
        }
        self = update_obj(self, args, force=True)
        for d in [self.bed_dir, self.bam_dir, self.bw_dir, self.matrix_dir]:
            if not os.path.exists(d):
                os.makedirs(d)

        
    # def split_bam_count(self, x):
    def count_bam_ss(self, x): # strand-specific
        ## split bam by strand: bam_count_dir
        xname = file_prefix(x)
        # (forward) -F 16
        c1_file = os.path.join(self.bam_dir, xname+'.count.fwd.txt')
        if os.path.exists(c1_file):
            with open(c1_file) as r:
                line = next(r).strip().split(',') # name,strand,count
                c1 = int(line[2])
        else:
            c1 = pysam.view('-c', '-F', '16', '-F', '4', '-@', '8', x) # fwd
            c1 = int(c1.strip())
            with open(c1_file, 'wt') as w:
                w.write(','.join([xname, 'fwd', str(c1)])+'\n')
        # (reverse) -f 16
        c2_file = os.path.join(self.bam_dir, xname+'.count.rev.txt')
        if os.path.exists(c2_file):
            with open(c2_file) as r:
                line = next(r).strip().split(',') # name,strand,count
                c2 = int(line[2])
        else:
            c2 = pysam.view('-c', '-f', '16', '-F', '4', '-@', '8', x) # fwd
            c2 = int(c2.strip())
            with open(c2_file, 'wt') as w:
                w.write(','.join([xname, 'fwd', str(c2)])+'\n')
        ## out
        return (c1, c2) # fwd, rev
    
    
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
    
    
    # strand-specific: False
    def bam2bw_ns(self, x, x_prefix=None): 
        args1 = self.__dict__.copy()
        if x_prefix is None:
            x_prefix = file_prefix(x)
        ## (forward)
        args = {
            'bam': x,
            'out_dir': self.bw_dir,
            'out_prefix': x_prefix,
            'normalizeUsing': 'RPGC',
        }
        args1.update(args)
        r = Bam2bw(**args1)
        r.run()
        return r.bw
    
        
    # strand-specific: True
    def bam2bw_ss(self, x, x_prefix=None):
        ## a. check norm-scale for fwd+rev # SE reads
        c1, c2 = self.count_bam_ss(x)
        args1 = self.__dict__.copy() 
        if x_prefix is None:
            x_prefix = file_prefix(x)
        # (forward)
        args_fwd = {
            'bam': x,
            'out_dir': self.bw_dir,
            'out_prefix': x_prefix+'_fwd',
            'filterRNAstrand': 'forward',
            'normalizeUsing': 'CPM', 
            'scaleFactor': c1/(c1+c2) # scaleFactor
        }
        args1.update(args_fwd)
        r1 = Bam2bw(**args1)
        r1.run()
        # (reverse)
        args2 = self.__dict__.copy() 
        args_rev = {
            'bam': x,
            'out_dir': self.bw_dir,
            'out_prefix': x_prefix+'_rev',
            'filterRNAstrand': 'reverse',
            'normalizeUsing': 'CPM', 
            'scaleFactor': c2/(c1+c2) # scaleFactor
        }
        args2.update(args_rev)
        r2 = Bam2bw(**args2)
        r2.run()
        ## output
        return (r1.bw, r2.bw)
        
        
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
        
        
    # to-do: multiple bam files
    # strand-specific: False
    def run_ns(self):
        """
        For non-stranded BAM files
        """
        # 1. bam to bw # to-do: in_parallel
        bw_list = [self.bam2bw_ns(i, j) for i,j in zip(self.bam_list, self.samplesLabel)]
        
        # 2. bw to matrix
        args2 = self.__dict__.copy()
        args2.update({
            'bw_list': bw_list,
            'out_dir': self.matrix_dir,
        })
        r2 = Bw2matrix(**args2)
        r2.run()
        
        # 3. copy matrix files
        src_dir_rel = os.path.relpath(os.path.dirname(r2.matrix), os.path.dirname(self.matrix))        
        src_rel = os.path.join(src_dir_rel, os.path.basename(r2.matrix))
        if not os.path.exists(self.matrix):
            os.symlink(src_rel, self.matrix)
        return self.matrix
        
        
    # strand-specific: True
    def run_ss(self):
        """
        For stranded BAM files
        """
        # out_dir_fwd, out_dir_rev = [os.path.join(self.out_dir, i) for i in ['fwd', 'rev']] # ?
        # 1. split regions file, by strand
        args1 = self.__dict__.copy()
        region_fwd_list = [self.split_bed(i, strand='+') for i in self.region_list]
        region_rev_list = [self.split_bed(i, strand='-') for i in self.region_list]
        
        # 2. bam to bw + matrix
        bw_list = [self.bam2bw_ss(i, j) for i,j in zip(self.bam_list, self.samplesLabel)]
        bw_fwd_list = [i[0] for i in bw_list]
        bw_rev_list = [i[1] for i in bw_list]
        # labels
        sl_fwd = [i+'_fwd' for i in self.samplesLabel]
        sl_rev = [i+'_rev' for i in self.samplesLabel]
        rl_fwd = [os.path.basename(i) for i in region_fwd_list] # 
        rl_rev = [os.path.basename(i) for i in region_rev_list] #
        
        # 3. bw to matrix
        args_bm = self.__dict__.copy()
        args_tmp = {
            'out_dir': self.matrix_dir,
            'out_prefix': self.out_prefix,
            'bw_list': bw_fwd_list + bw_rev_list,
            'region_list': region_fwd_list + region_rev_list,
            'samplesLabel': sl_fwd + sl_rev,
            'regionsLabel': rl_fwd + rl_rev,
        }
        args_bm.update(args_tmp)
        s1 = Bw2matrix(**args_bm)
        s1.run()
        
        # 4. subset matrix
        ## sense
        sens1 = os.path.join(self.matrix_dir, self.out_prefix+'.sens_1.mat.gz')
        sens2 = os.path.join(self.matrix_dir, self.out_prefix+'.sens_2.mat.gz')
        self.subset_matrix(s1.matrix, sens1, samples=sl_fwd, groups=rl_fwd) # fwd+fwd
        self.subset_matrix(s1.matrix, sens2, samples=sl_rev, groups=rl_rev) # rev+rev
        ## antisense
        anti1 = os.path.join(self.matrix_dir, self.out_prefix+'.anti_1.mat.gz')
        anti2 = os.path.join(self.matrix_dir, self.out_prefix+'.anti_2.mat.gz')
        self.subset_matrix(s1.matrix, anti1, samples=sl_fwd, groups=rl_rev) # fwd+rev
        self.subset_matrix(s1.matrix, anti2, samples=sl_rev, groups=rl_fwd) # rev+fwd
        ## merge
        # self.rbind_matrix(x=[sens1, sens2], o=self.matrix_sens)
        # self.rbind_matrix(x=[anti1, anti2], o=self.matrix_anti)
        Matrix_rbind(m=[sens1, sens2], o=self.matrix_sens).run()
        Matrix_rbind(m=[anti1, anti2], o=self.matrix_anti).run()
        return (self.matrix_sens, self.matrix_anti)
        
        
    def run(self):
        print(self.check_files())
        if self.strand_specific:
            self.run_ss()
        else:
            self.run_ns()


def get_args():
    example = ' '.join([
        '$ python bam2matrix.py', 
        '-b f1.bam f2.bam -r g1.bed g2.bed -o out_dir --out-prefix metaplot', 
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
        help='bw files')
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
    Bam2matrix(**args).run()

    
if __name__ == '__main__':
    main()

# EOF