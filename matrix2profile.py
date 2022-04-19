#!/usr/bin/env python

"""
Generate profile plot 

$ plotProfile \
  -m input.mat -o out.png \
  --samplesLabel A B C --regionsLabel gene1 gene2 \
  --colors black lightblue --yMin 0 --yMax 0.4 \
  --perGroup
  # --numPlotsPerRow 2
  # --plotTitle
  # --plotHeight 4
  # --plotWidth 1
  # --clusterUsingSamples 1
  # --refPointLabel TSS
  # --startLabel TSS --endLabel TES
"""

import os
import sys
import pathlib
import argparse
import shutil
from matplotlib import colors

from utils import (
    make_config, update_obj, load_yaml, dump_yaml, file_abspath, file_prefix,
    fix_label, fix_bw, is_valid_bam, is_valid_bigwig, is_valid_bed, log,
    load_matrix_header
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
      # --numPlotsPerRow 2
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
        self.init_args()
        self.init_files()
        self.cmd = self.get_cmd()
        # save config
        dump_yaml(self.__dict__, self.config)


    def init_args(self):
        if not hasattr(self, 'colors'):
            self.colors = None
        # check input files + labels
        self.matrix = file_abspath(self.matrix)
        d = load_matrix_header(self.matrix)
        if d is None:
            raise ValueError('unable to read matrix: {}'.format(self.matrix))
        # check if reference-point or scale-regions
        bd = d.get('body', [0])
        rf = d.get('ref point', [None])
        self.is_refpoint = isinstance(rf[0], str) and bd[0] == 0
        # update labels
        d_samples_label = d.get('sample_labels', [None])
        d_regions_label = d.get('group_labels', [None])
        if isinstance(self.samplesLabel, list):
            self.is_valid_sl = len(d_samples_label) == len(self.samplesLabel)
        else:
            self.is_valid_sl = False
        if isinstance(self.regionsLabel, list):
            self.is_valid_rl = len(d_regions_label) == len(self.regionsLabel)
        else:
            self.is_valid_rl = False


    def init_files(self):
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
        # basic
        args_basic = ' '.join([
            '{}'.format(shutil.which('plotProfile')),
            '--matrixFile {}'.format(self.matrix),
            '--outFileName {}'.format(self.profile_file),
            '--dpi {}'.format(self.dpi),
        ])
        if self.is_valid_sl:
            args_basic += ' --samplesLabel {}'.format(' '.join(self.samplesLabel))
        if self.is_valid_rl:
            args_basic += ' --regionsLabel {}'.format(' '.join(self.regionsLabel))
        # type
        if self.is_refpoint:
            args_type = '--refPointLabel {}'.format(self.refPointLabel)
        else:
            args_type = ' '.join([
                '--startLabel {}'.format(self.startLabel),
                '--endLabel {}'.format(self.endLabel),
            ])
        # extra
        args_extra = '--perGroup' if self.perGroup else ''
        if isinstance(self.colors, str):
            cc = self.colors
        elif isinstance(self.colors, list):
            cc = ' '.join(self.colors)
        else:
            cc = None
        if cc is not None:
            args_extra += ' --colors {}'.format(cc)
        if isinstance(self.yMin, float):
            args_extra += ' --yMin {}'.format(self.yMin)
        if isinstance(self.yMax, float):
            args_extra += ' --yMax {}'.format(self.yMax)
        # log
        args_extra += ' 1> {}'.format(self.stdout)
        args_extra += ' 2> {}'.format(self.stderr)
        # cmd
        return ' '.join([args_basic, args_type, args_extra])


    def run(self):
        # save command
        with open(self.profile_cmd, 'wt') as w:
            w.write(self.cmd+'\n')
        # run
        print(self.overwrite)
        if os.path.exists(self.profile_file) and not self.overwrite:
            # if re-cal required, remove the old file
            log.info('plotProfile() skipped, file exists: {}'.format(self.profile_file))
        else:
            log.info('run plotProfile: {}'.format(self.cmd))
            os.system(self.cmd)


def get_args():
    example = ' '.join([
        '$ plotProfile',
        '-m input.mat -o out.png',
        '--samplesLabel A B C --regionsLabel gene1 gene2',
        '--colors black lightblue --yMin 0 --yMax 0.4',
        '--perGroup',
    ])
    parser = argparse.ArgumentParser(prog='bw2matrix.py',
                                     description='bw2matrix',
                                     epilog=example,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-m', dest='matrix', required=True,
        help='matrix file, by computeMatrix ')
    parser.add_argument('-o', dest='out', required=True,
        help='filename of output')    
    parser.add_argument('--averageType', default='mean',
        choices=['mean', 'median', 'max', 'min', 'sum', 'region_length'],
        help='Which method should be used for sorting, default: [mean]')
    parser.add_argument('--plotType', default='lines',
        choices=['lines', 'fill', 'se', 'std', 'overlapped_lines', 'heatmap'],
        help='type of the plot, default: [lines]')
    parser.add_argument('--colors', nargs='+', default=None,
        help='colors for the lines, default: [None] auto')
    parser.add_argument('--regionsLabel', nargs='+', default=None,
        help='labels for regions in plot, defautl: [None] auto')
    parser.add_argument('--samplesLabel', nargs='+', default=None,
        help='labels for samples in plot, default: [None] auto')
    parser.add_argument('--startLabel', default='TSS',
        help='start label, default: [TSS]')
    parser.add_argument('--endLabel', default='TES',
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
    args.update({
        'out_dir': os.path.dirname(args['out']),
        'out_prefix': file_prefix(args['out']),
    })
    Matrix2profile(**args).run()

    
if __name__ == '__main__':
    main()

# EOF