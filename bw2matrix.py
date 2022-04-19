#!/usr/bin/env python

"""
Convert bigWig to matrix

1. general purpose
- sample_list
- region_list
- samplesLabel
- matrix # output

## scale-regions
- b 2000
- a 2000
- m 2000
- binSize 100
- p 8
- startLabel TSS
- endLabel TES

## reference-point
- b 2000
- a 2000
- binSize 100
- p 8
- referencePoint  TSS

$ computeMatrix scale-regions \
  -R in.bed -S a.bw -o out.mat.gz \
  -b 2000 -a 2000 -m 2000 --binSize 100 -p 8 \
  --unscaled5prime 0 --unscaled3prime 0 --skipZeros \
  --averageTypeBins mean --blackListFileName bl.bed \
  --outFileSortedRegions sorted.bed \
  --startLabel TSS --endLabel TES
  # --missingDataAsZero
  # --samplesLabel
  # --smartLabels

$ computeMatrix reference-point \
  -R in.bed -S a.bw -o out.mat.gz \
  -b 2000 -a 2000 --binSize 100 -p 8 \
  --unscaled5prime 0 --unscaled3prime 0 --skipZeros \
  --averageTypeBins mean --blackListFileName bl.bed \
  --outFileSortedRegions sorted.bed
"""

import os
import sys
import pathlib
import argparse
import shutil

from utils import (
    make_config, update_obj, dump_yaml, file_abspath, file_prefix, symlink_file,
    fix_label, fix_bw, is_valid_bigwig, is_valid_bed, is_valid_file, 
    load_matrix_header, log
)
from matrix_rbind import Matrix_rbind

################################################################################
## 2. BigWig to matrix
class Bw2matrix(object):
    """
    Strand-specific: True/False
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
            self.bw_list = None
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
            self.bw_fwd_list = self.bw_rev_list = None
            self.samplesLabel = fix_label(self.bw_rwd_list, self.samplesLabel)
        else:
            raise ValueError('bw_list, bw_fwd_list, bw_rev_list required:')
        if not isinstance(self.out_prefix, str):
            self.out_prefix = 'metaplot'
            log.warning('out_prefix failed, use "metaplot" instead')
        # type
        if not self.matrix_type in ['scale-regions', 'reference-point']:
            self.matrix_type = 'scale-regions'
            

    def init_files(self):
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.out_prefix)
        if not os.path.exists(self.project_dir):
            os.makedirs(self.project_dir)


    def run(self):
        fun = Bw2matrix_ss if self.strand_specific else Bw2matrix_ns
        args = self.__dict__.copy()
        r = fun(**args)
        return r.run()


class Bw2matrix_ns(object):
    """
    Strand-specific: False
    """
    def __init__(self, **kwargs):
        c = make_config(**kwargs)
        self = update_obj(self, c, force=True)
        self.init_args()
        self.init_files()
        self.cmd = self.get_cmd()
        # save config
        dump_yaml(self.__dict__, self.config)


    def init_args(self):
        pass # in Bw2matrix().init_args()
#         # region files (list)
#         if is_valid_file(self.region_list, is_valid_bed):
#             if isinstance(self.region_list, str):
#                 self.region_list = [self.region_list]
#             self.region_list = list(map(file_abspath, self.region_list))
#             self.regionsLabel = fix_label(self.region_list, self.regionsLabel)
#         else:
#             raise ValueError('illegal, region_list={}'.format(
#                 self.region_list
#             ))
#         # bw files
#         if is_valid_file(self.bw_list, is_valid_bigwig):
#             self.strand_specific = False
#             self.samplesLabel = fix_label(self.bw_list, self.samplesLabel)
#         else:
#             raise ValueError('illegal, bw_list={}'.format(
#                 self.bw_list
#             ))
#         if not isinstance(self.out_prefix, str):
#             self.out_prefix = 'metaplot'
#             log.warning('out_prefix failed, use "metaplot" instead')
#         # type
#         if not self.matrix_type in ['scale-regions', 'reference-point']:
#             self.matrix_type = 'scale-regions'


    def init_files(self):
#         if not isinstance(self.out_dir, str):
#             self.out_dir = str(pathlib.Path.cwd())
#         self.out_dir = file_abspath(self.out_dir)
#         self.project_dir = os.path.join(self.out_dir, self.out_prefix)
        prefix = os.path.join(self.project_dir, self.out_prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.yaml'),
            'matrix': prefix+'.mat.gz',
            'matrix_value': prefix+'.mat.tab',
            'sorted_regions_file': prefix+'.sortedRegions.bed',
            'stdout': prefix+'.computeMatrix.stdout',
            'stderr': prefix+'.computeMatrix.stderr',
            'matrix_cmd': os.path.join(self.project_dir, prefix+'.computeMatrix.sh'),
        }
#         if not os.path.exists(self.project_dir):
#             os.makedirs(self.project_dir)
        self = update_obj(self, args, force=True)


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


    def get_cmd(self):
        args_basic = ' '.join([
            '{}'.format(shutil.which('computeMatrix')),
            '{}'.format(self.matrix_type)
        ])
        # I/O
        args_io = ' '.join([
            '-S {}'.format(' '.join(self.bw_list)),
            '-R {}'.format(' '.join(self.region_list)),
            '-o {}'.format(self.matrix),
            '--outFileNameMatrix {}'.format(self.matrix_value),
            '--outFileSortedRegions {}'.format(self.sorted_regions_file),
        ])
        # size
        args_size = ' '.join([
            '--skipZeros',
            '--missingDataAsZero',
            '--numberOfProcessors {}'.format(self.numberOfProcessors),
            '--binSize {}'.format(self.binSize),
            '--beforeRegionStartLength {}'.format(self.beforeRegionStartLength),
            '--afterRegionStartLength {}'.format(self.afterRegionStartLength),
            '--sortRegions {}'.format(self.sortRegions),
            '--sortUsing {}'.format(self.sortUsing),
        ])
        if isinstance(self.sortUsingSamples, str):
            args_size += ' --sortUsingSamples {}'.format(self.sortUsingSamples)
        # type
        if self.matrix_type == 'reference-point':
            args_type = ' '.join([
                '--referencePoint {}'.format(self.referencePoint),
            ])
        else:
            args_type = ' '.join([
                '--regionBodyLength {}'.format(self.regionBodyLength),
                '--startLabel {}'.format(self.startLabel),
                '--endLabel {}'.format(self.endLabel),
                '--unscaled5prime {}'.format(self.unscaled5prime),
                '--unscaled3prime {}'.format(self.unscaled3prime),
            ])
        # extra
        args_extra = ' '.join([
            ' --samplesLabel {}'.format(' '.join(self.samplesLabel)),
            '1> {}'.format(self.stdout),
            '2> {}'.format(self.stderr),
        ])
        # out
        return ' '.join([args_basic, args_io, args_size, args_type, args_extra])


    def run(self):
        print(self.check_files())
        # save command
        with open(self.matrix_cmd, 'wt') as w:
            w.write(self.cmd+'\n')
        # run
        if os.path.exists(self.matrix) and not self.overwrite:
            log.info('Bw2matrix() skipped, file exists: {}'.format(self.matrix))
        else:
            log.info('run computeMatrix: {}'.format(self.matrix))
            os.system(self.cmd)
        # check output
        if not os.path.exists(self.matrix):
            log.error('Bw2matrix() failed: {}'.format(self.stderr))
        return self.matrix


class Bw2matrix_ss(object):
    """
    Strand-specific: True
    """
    def __init__(self, **kwargs):
        c = make_config(**kwargs)
        self = update_obj(self, c, force=True)
        self.init_args()
        self.init_files()


    def init_args(self):
        pass # in Bw2matr
#         # bw files
#         if is_valid_file(self.bw_fwd_list, is_valid_bigwig) and \
#             is_valid_file(self.bw_rev_list, is_valid_bigwig):
#             self.strand_specific = True
#             tag = 0
#             for f,r in zip(self.bw_fwd_list, self.bw_rev_list):
#                 f1 = file_prefix(f).strip('_fwd')
#                 r1 = file_prefix(r).strip('_rev')
#                 if not f1 == r1:
#                     tag += 1
#                     print('bw fwd/rev not matched: {} {}'.format(f, r))
#             if tag > 0:
#                 raise ValueError('bw fwd/rev files illegal')
#             # update samplesLabel
#             self.samplesLabel = fix_label(self.bw_fwd_list, self.samplesLabel)
#             self.samplesLabel = [i.rstrip('_fwd') for i in self.samplesLabel]
#         else:
#             raise ValueError('bw_list, bw_fwd_list, bw_rev_list required:')
#         if not isinstance(self.out_prefix, str):
#             self.out_prefix = 'metaplot'
#             log.error('out_prefix failed, use "metaplot" instead')


    def init_files(self):
#         if not isinstance(self.out_dir, str):
#             self.out_dir = str(pathlib.Path.cwd())
#         self.out_dir = file_abspath(self.out_dir)
#         self.project_dir = os.path.join(self.out_dir, self.out_prefix)
        prefix = os.path.join(self.project_dir, self.out_prefix)
        args = {
#             'config': os.path.join(self.project_dir, 'config.yaml'),
            'matrix_sens': prefix+'_sens.mat.gz',
            'matrix_anti': prefix+'_anti.mat.gz',
#             'matrix_value_sens': prefix+'_sens.mat.tab',
#             'matrix_value_anti': prefix+'_anti.mat.tab',
#             'sorted_regions_file': prefix+'.sortedRegions.bed',
#             'stdout': prefix+'.computeMatrix.stdout',
#             'stderr': prefix+'.computeMatrix.stderr',
#             'matrix_cmd_sens': os.path.join(self.project_dir, prefix+'_sens.computeMatrix.sh'),
#             'matrix_cmd_anti': os.path.join(self.project_dir, prefix+'_anti.computeMatrix.sh'),
        }
        if not os.path.exists(self.project_dir):
            os.makedirs(self.project_dir)
        self = update_obj(self, args, force=True)


    def split_bed(self, x, strand='+'):
        """
        Split BED file by strand
        * : original
        + : _fwd
        - : _rev
        """
        xname = file_prefix(x)
        t = '_rev' if strand == '-' else '_fwd' if strand == '+' else ''
        out = os.path.join(self.project_dir, xname+t+'.bed')
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
            symlink_file(x, o) # copy
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
                try:
                    os.system(cmd) # try-except ?
                except:
                    log.error('subset_matrix() failed, check: {}'.format(cmd))
        # check output
        if not os.path.exists(o):
            log.error('subset_matrix() failed, file not exists: {}'.format(o))


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
            # 'out_dir': self.out_dir,
            # 'out_prefix': self.out_prefix,
            'bw_list': self.bw_fwd_list + self.bw_rev_list,
            'region_list': region_fwd_list + region_rev_list,
            'samplesLabel': sl_fwd + sl_rev,
            'regionsLabel': rl_fwd + rl_rev,
        })
        r2 = Bw2matrix_ns(**args2)
        r2.run()

        # 3. subset matrix
        ## sense
        sens1 = os.path.join(self.project_dir, self.out_prefix+'_sens_1.mat.gz')
        sens2 = os.path.join(self.project_dir, self.out_prefix+'_sens_2.mat.gz')
        self.subset_matrix(r2.matrix, sens1, samples=sl_fwd, groups=rl_fwd) # fwd+fwd
        self.subset_matrix(r2.matrix, sens2, samples=sl_rev, groups=rl_rev) # rev+rev
        ## antisense
        anti1 = os.path.join(self.project_dir, self.out_prefix+'_anti_1.mat.gz')
        anti2 = os.path.join(self.project_dir, self.out_prefix+'_anti_2.mat.gz')
        self.subset_matrix(r2.matrix, anti1, samples=sl_fwd, groups=rl_rev) # fwd+rev
        self.subset_matrix(r2.matrix, anti2, samples=sl_rev, groups=rl_fwd) # rev+fwd
        # 4. merge
        Matrix_rbind(m=[sens1, sens2], o=self.matrix_sens).run()
        Matrix_rbind(m=[anti1, anti2], o=self.matrix_anti).run()
        return (self.matrix_sens, self.matrix_anti)


    def run(self):
        return self.bw2matrix_ss()


def get_args():
    example = ' '.join([
        '$ computeMatrix scale-regions',
        '-R in.bed -S a.bw -o out.mat.gz',
        '-b 2000 -a 2000 -m 2000 --binSize 100 -p 8',
        '--unscaled5prime 0 --unscaled3prime 0 --skipZeros',
        '--averageTypeBins mean --blackListFileName bl.bed',
        '--outFileSortedRegions sorted.bed',
        '--startLabel TSS --endLabel TES',
    ])
    parser = argparse.ArgumentParser(prog='bw2matrix.py',
                                     description='bw2matrix',
                                     epilog=example,
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
    parser.add_argument('--out-prefix', dest='out_prefix', default='bw2matrix',
        help='prefix for output files, default: [bw2matrix]')
    parser.add_argument('-t', '--matrix-type', dest='matrix_type',
        default='scale-regions', choices=['scales-regions', 'reference-point'],
        help='choose the matrix type, default: [scale-regions]')
    parser.add_argument('-sl', '--samplesLabel', default=None,
        help='samples label')
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
    Bw2matrix(**args).run()


if __name__ == '__main__':
    main()

# EOF