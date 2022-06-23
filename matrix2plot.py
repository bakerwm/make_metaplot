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
import argparse
from utils import (
    make_config, update_obj, log,
)
from matrix2profile import Matrix2profile
from matrix2heatmap import Matrix2heatmap


class Matrix2plot(object):
    """
    Generate plots using matrix
    """
    def __init__(self, **kwargs):
        c = make_config(**kwargs)
        self = update_obj(self, c, force=True)


    def run_profile(self):
        args = self.__dict__.copy()
        args.update({
            'out_dir': os.path.join(self.out_dir, '3.matrix2profile')
        })
        return Matrix2profile(**args).run()


    def run_heatmap(self):
        args = self.__dict__.copy()
        args.update({
            'out_dir': os.path.join(self.out_dir, '4.matrix2heatmap')
        })
        return Matrix2heatmap(**args).run()

    
    def run(self):
        self.run_profile()
        self.run_heatmap()


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