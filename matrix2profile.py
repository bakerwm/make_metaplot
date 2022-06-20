#!/usr/bin/env python
#-*- encoding:utf-8 -*-
"""
Generate profile plot 

$ plotProfile \
  -m input.mat -o out.png \
  --samplesLabel A B C --regionsLabel gene1 gene2 \
  --colors black lightblue --yMin 0 --yMax 0.4 \
  --perGroup
  --numPlotsPerRow 2
  # --plotTitle
  # --plotHeight 4
  # --plotWidth 1
  # --clusterUsingSamples 1
  # --refPointLabel TSS
  # --startLabel TSS --endLabel TES
"""

import os
import pathlib
import argparse
import shutil
from matplotlib import colors

from utils import (
    make_config, update_obj, dump_yaml, file_abspath, file_prefix,
    fix_label, is_valid_bed, is_valid_file, log,
    load_matrix, load_matrix_header
)


## profile, heatmap, ...
class Matrix2profile(object):
    """
    Example:
    $ plotProfile \
      -m input.mat -o out.png \
      --samplesLabel A B C --regionsLabel gene1 gene2 \
      --colors black lightblue --yMin 0 --yMax 0.4 \
      --perGroup
      --numPlotsPerRow 2
      # --plotTitle
      # --plotHeight 4
      # --plotWidth 1
      # --clusterUsingSamples 1
      # --refPointLabel TSS
      # --startLabel TSS --endLabel TES
    """
    def __init__(self, **kwargs):
        c = make_config(**kwargs)
        self = update_obj(self, c, force=True)
        self.update_args() # set default to None
        self.update_labels()
        self.update_colors()
        self.init_files()
        self.cmd = self.get_cmd()
        dump_yaml(self.__dict__, self.config) # save config


    def basic_args(self):
        """
        default arguments [38]
        to-do: outFileSortedRegions, outFileNameData
        """        
        alist = [
            'startLabel', 'endLabel', 'refPointLabel', 'samplesLabel', 
            'regionsLabel', 'plotTitle', 'yAxisLabel', 'labelRotation',
            'colors', 'numPlotsPerRow', 'clusterUsingSamples', 
            'plotHeight', 'plotWidth', 'plotType',
            'yMin', 'yMax', 'legendLocation',
            'kmeans', 'hclust', 'silhouette', 'dpi', 'averageType',
            'outFileSortedRegions', 'outFileNameData',
        ]
        return alist


    def update_args(self):
        """
        assign None to arguments, if not exists
        """
        alist = self.basic_args()
        d = {i:getattr(self, i, None) for i in alist}
        self = update_obj(self, d, force=True) # update        
        self.whatToShow = '"{}"'.format(self.whatToShow) # update whatToShow


    def update_labels(self):
        """
        convert to str, or keep None
        1. samplesLabel, regionsLabel
        2. startLabel, endLabel / refPointLabel 
        3. yAxisLabel
        4. plotTitle, legendLocation
        """
        # load matrix
        mh = load_matrix(self.matrix, header_only=True)
        mh_sl = mh.get('sample_labels', [None])
        mh_rl = mh.get('group_labels', [None])
        is_refPoint = mh.get('body', 0) == 0 # body == 0
        # 1. samplesLabel, regionsLabel
        if isinstance(self.samplesLabel, list):
            k1 = len(self.samplesLabel) == len(mh_sl)
            # raise error
            self.samplesLabel = ' '.join(self.samplesLabel)
        if isinstance(self.regionsLabel, list):
            k2 = len(self.regionsLabel) == len(mh_rl)
            # raise error
            self.regionsLabel = ' '.join(self.regionsLabel)
        # 2. startLabel, endLabel or refPointLabel
        if is_refPoint:
            self.startLabel, self.endLabel = [None, None]
        else:
            self.refPointLabel = None
        # 3. yAxisLabel
        # 4. plotTitle, legendLocation


    def update_colors(self):
        """
        colorMap: Reds Blues
        colorList: white,red white,yellow,blue
        colorNumber: int
        """
        if isinstance(self.colors, list):
            self.colors = ' '.join(self.colors)


    def init_files(self):
        self.matrix = file_abspath(self.matrix)
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.out_prefix)
        prefix = os.path.join(self.project_dir, self.out_prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.yaml'),
            'profile_cmd': os.path.join(self.project_dir, prefix+'.plotProfile.sh'),
            'profile_file': prefix+'.plotProfile.pdf',
            'stdout': prefix+'.plotProfile.stdout',
            'stderr': prefix+'.plotProfile.stderr',
        }
        if not os.path.exists(self.project_dir):
            os.makedirs(self.project_dir)
        self = update_obj(self, args, force=True)


    def get_cmd(self):
        """
        construct arguments to command line
        """
        alist = self.basic_args()
        # args = self.__dict__.copy() #
        args = {i:getattr(self, i, None) for i in alist}
        dlist = ['--{} {}'.format(k, v) for k,v in args.items() if v is not None]
        dline = ' '.join(dlist) # to cmd line
        if self.perGroup:
            dline += ' --perGroup'
        # main args
        cmd = ' '.join([
            '{}'.format(shutil.which('plotProfile')),
            '--matrixFile {}'.format(self.matrix),
            '--outFileName {}'.format(self.profile_file),
            dline,            
            '1> {}'.format(self.stdout),
            '2> {}'.format(self.stderr),
        ])
        return cmd


    def run(self):
        # save command
        with open(self.profile_cmd, 'wt') as w:
            w.write(self.cmd+'\n')
        # run
        if os.path.exists(self.profile_file) and not self.overwrite:
            # if re-cal required, remove the old file
            log.info('plotProfile() skipped, file exists: {}'.format(self.profile_file))
        else:
            log.info('run plotProfile: {}'.format(self.cmd))
            os.system(self.cmd)
        # check output
        if not os.path.exists(self.profile_file):
            log.error('plotProfile() failed, file not found: {}'.format(self.profile_file))


def get_args():
    example = ' '.join([
        'Example: \n',
        '$ python matrix2profile.py',
        '-m aaa/bw2matrix/bw2matrix_sens.mat.gz -o aaa',
        '-op metaplot --plotType se ',
        '--colors black red -p 2 -rl genes',
    ])
    parser = argparse.ArgumentParser(
        prog='bw2matrix.py', description='bw2matrix', epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-m', dest='matrix', required=True,
        help='matrix file, by computeMatrix ')
    parser.add_argument('-o', dest='out_dir', required=False,
        help='directory to save bigWig file')
    parser.add_argument('-op', '--out-prefix', dest='out_prefix', default='metaplot',
        help='prefix for output files, default: [metaplot]')
    parser.add_argument('--plotType', default='lines',
        choices=['lines', 'fill', 'se', 'std', 'overlapped_lines', 'heatmap'],
        help='type of the plot, default: [lines]')
    parser.add_argument('--colors', nargs='+', default=None,
        help='colors for the lines, default: [None] auto')
    parser.add_argument('--averageType', default='mean',
        choices=['mean', 'median', 'max', 'min', 'sum', 'region_length'],
        help='Which method should be used for sorting, default: [mean]')
    parser.add_argument('-rl', '--regionsLabel', nargs='+', default=None,
        help='labels for regions in plot, defautl: [None] auto')
    parser.add_argument('-sl', '--samplesLabel', nargs='+', default=None,
        help='labels for samples in plot, default: [None] auto')
    parser.add_argument('-st', '--startLabel', default='TSS',
        help='start label, default: [TSS]')
    parser.add_argument('-ed', '--endLabel', default='TES',
        help='end label, default: [TES]')
    parser.add_argument('--refPointLabel', default='TSS',
        help='refPointLabel label, default: [TSS]')
    parser.add_argument('--yMin', type=float, default=None,
        help='Minimum value for the Y-axis')
    parser.add_argument('--yMax', type=float, default=None,
        help='Maximum value for the Y-axis')
    parser.add_argument('--perGroup', action='store_true',
        help='plot all samples by group')
    parser.add_argument('-p', dest='numberOfProcessors', type=int, default=4, 
        help='number of processors, default: [4]')
    parser.add_argument('-O', '--overwrite', dest='overwrite', action='store_true',
        help='Overwrite output file')
    return parser


def main():
    args = vars(get_args().parse_args())
    Matrix2profile(**args).run()

    
if __name__ == '__main__':
    main()

# EOF