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
    make_config, update_obj, load_yaml, dump_yaml, file_abspath, file_prefix,
    fix_label, fix_bw, is_valid_bigwig, is_valid_bed, is_valid_file, log
)

################################################################################
## 2. BigWig to matrix
class Bw2matrix(object):
    """
    Example:
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
    def __init__(self, **kwargs):
        c = make_config(**kwargs)
        self = update_obj(self, c, force=True)
        self.init_args()
        self.init_files()
        self.cmd = self.get_cmd()
        # save config
        dump_yaml(self.__dict__, self.config)


    def init_args(self):
        # check input files + labels
        if isinstance(self.bw_list, str):
            self.bw_list = [self.bw_list]
        if isinstance(self.region_list, str):
            self.region_list = [self.region_list]
        # if isinstance(self.bw_list, list):
        self.bw_list = list(map(file_abspath, self.bw_list))
        self.samplesLabel = fix_label(self.bw_list, self.samplesLabel)
        # if isinstance(self.region_list, list):
        self.region_list = list(map(file_abspath, self.region_list))
        self.regionsLabel = fix_label(self.region_list, self.regionsLabel)
        # type
        if not self.matrix_type in ['scale-regions', 'reference-point']:
            self.matrix_type = 'scale-regions'
        # output
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)


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
        # valid files
        c1 = is_valid_file(self.bw_list, is_valid_bigwig)
        c2 = is_valid_file(self.region_list, is_valid_bed)
        if not c1 and c2:
            msg2 = '\n'.join([
                'Check the following files:',
                '-'*80,
                '>> score files:',
                '\n'.join([
                '{} : {}'.format(is_valid_bigwig(i), i) for i in self.bw_list
                ]),
                '-'*30,
                '>> region files:',
                '\n'.join([
                '{} : {}'.format(is_valid_bed(i), i) for i in self.region_list
                ]),
                '-'*80,
            ])
            raise Exception(msg2)
        return msg


    def init_files(self):
        self.project_dir = os.path.join(self.out_dir, self.out_prefix)
        prefix = os.path.join(self.project_dir, self.out_prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.yaml'),
            'matrix': prefix+'.mat.gz',
            'matrix_value': prefix+'.mat.tab',
            'sorted_regions_file': prefix+'.sortedRegions.bed',
            'stdout': prefix+'.computeMatrix.stdout',
            'stderr': prefix+'.computeMatrix.stderr',
            'matrix_cmd': os.path.join(self.project_dir, prefix+'.computeMatrix.sh')
        }
        if not os.path.exists(self.project_dir):
            os.makedirs(self.project_dir)
        self = update_obj(self, args, force=True)


    def get_cmd(self):
        # basic
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
        args_extra = ''
        if isinstance(self.samplesLabel, list):
            args_extra += ' --samplesLabel {}'.format(' '.join(self.samplesLabel))
#         if isinstance(self.regionsLabel, list):
#             args_extra += ' --regionsLabel {}'.format(' '.join(self.regionsLabel))
        # log
        args_extra += ' 1> {}'.format(self.stdout)
        args_extra += ' 2> {}'.format(self.stderr)
        # cmd
        return ' '.join([args_basic, args_io, args_size, args_type, args_extra])


    def run(self):
        print(self.check_files())
        # save command
        with open(self.matrix_cmd, 'wt') as w:
            w.write(self.cmd+'\n')
        # run
        if os.path.exists(self.matrix) and not self.overwrite:
            # if re-cal required, remove the old file
            log.info('Bw2matrix() skipped, file exists: {}'.format(self.matrix))
        else:
            log.info('run computeMatrix: {}'.format(self.matrix))
            os.system(self.cmd)
        # check output
        if not os.path.exists(self.matrix):
            log.error('Bw2matrix() failed: {}'.format(self.stderr))
            

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
    parser.add_argument('-b', dest='bw_list', nargs='+', required=True,
        help='bw files')
    parser.add_argument('-r', dest='region_list', nargs='+', required=True,
        help='region files')
    parser.add_argument('-o', dest='out_dir', required=True,
        help='directory to save bigWig file')
    parser.add_argument('--out-prefix', dest='out_prefix', default='bw2matrix',
        help='prefix for output files, default: [bw2matrix]')
    parser.add_argument('-t', '--matrix-type', dest='matrix_type', 
        default='scale-regions', choices=['scales-regions', 'reference-point'],
        help='choose the matrix type, default: [scale-regions]')    
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
    Bw2matrix(**args).run()

    
if __name__ == '__main__':
    main()

# EOF