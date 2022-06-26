#!/usr/bin/env python
#-*- encoding:utf-8 -*-
"""
Generate heatmap
"""

import os
import pathlib
import argparse
import shutil
from matplotlib import colors
from utils import (
    make_config, update_obj, dump_yaml, file_prefix, file_abspath, log, 
    load_matrix, fix_out_dir,
)
from parse_args import  add_plot_parser, get_plot_args


class Matrix2heatmap(object):
    """
    Example:
    $ plotHeatmap \
      -m input.mat -o out.png \
      --whatToShow "plot, heatmap and colorbar" \
      --heatmapHeight 15 --heatmapWidth 4 \
      --colorList red,white,blue \
      --samplesLabel WT 15m 30m 120m --regionsLabel PAS-both
    """
    def __init__(self, **kwargs):
        # c = make_config(**kwargs)
        c = get_plot_args(**kwargs)
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
        to-do: outFileSortedRegions, outFileNameMatrix
        """        
        alist = [
            'startLabel', 'endLabel', 'refPointLabel', 'samplesLabel', 
            'regionsLabel', 'plotTitle', 'xAxisLabel', 'yAxisLabel',
            'legendLocation', 'colorMap', 'colorList', 'colorNumber',
            'outFileSortedRegions', 'outFileNameMatrix',
            'plotType', 'sortRegions', 'sortUsing', 'averageTypeSummaryPlot', 
            'alpha', 'zMin', 'zMax', 'yMin', 'yMax', 'heatmapHeight', 
            'heatmapWidth', 'whatToShow', 'labelRotation', 
            'kmeans', 'hclust', 'silhouette', 'sortUsingSamples',
            'clusterUsingSamples', 'missingDataColor',
            'boxAroundHeatmaps', 'interpolationMethod', 'dpi'
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
        if not isinstance(self.prefix, str):
            self.prefix = file_prefix(self.matrix)


    def update_labels(self):
        """
        convert to str, or keep None
        1. samplesLabel, regionsLabel
        2. startLabel, endLabel / refPointLabel 
        3. xAxisLabel, yAxisLabel
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
        # 3. xAxisLabel, yAxisLabel
        if isinstance(self.xAxisLabel, str):
            self.xAxisLabel = '"{}"'.format(self.xAxisLabel)
        if isinstance(self.yAxisLabel, str):
            self.yAxisLabel = '"{}"'.format(self.yAxisLabel)
        # 4. plotTitle, legendLocation


    def update_colors(self):
        """
        colors: 
        colorMap: Reds Blues
        colorList: white,red white,yellow,blue
        colorNumber: int
        """
        if isinstance(self.colorMap, list):
            self.colorMap = ' '.join(self.colorMap)
        if isinstance(self.colorList, list):
            self.colorList = ' '.join(self.colorList)
        if not isinstance(self.colorNumber, int):
            self.colorNumber = None
        # update colors: add " " to color names
        if isinstance(self.colors, list):
            self.colors = ' '.join(self.colors)
        if isinstance(self.colors, str):
            cc = ['"{}"'.format(i) for i in self.colors.split()]
            self.colors = ' '.join(cc)


    def init_files(self):
        self.matrix = file_abspath(self.matrix)
        self.out_dir = fix_out_dir(self.out_dir)
        self.project_dir = self.out_dir
        # if not isinstance(self.out_dir, str):
        #     self.out_dir = str(pathlib.Path.cwd())
        # self.out_dir = file_abspath(self.out_dir)
        prefix = os.path.join(self.project_dir, self.prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.yaml'),
            'heatmap_cmd': os.path.join(self.project_dir, prefix+'.plotHeatmap.sh'),
            'heatmap_file': prefix+'.plotHeatmap.pdf',
            'stdout': prefix+'.plotHeatmap.stdout',
            'stderr': prefix+'.plotHeatmap.stderr',
        }
        if not os.path.exists(self.project_dir):
            os.makedirs(self.project_dir)
        self = update_obj(self, args, force=True)


    def get_cmd(self):
        """
        construct arguments to command line:
        ['perGroup', 'linesAtTickMarks']
        """
        alist = self.basic_args()
        # args = self.__dict__.copy() #
        args = {i:getattr(self, i, None) for i in alist}
        dlist = ['--{} {}'.format(k, v) for k,v in args.items() if v is not None]
        bb = ['perGroup', 'linesAtTickMarks'] # no arguments
        bba = ['--'+i for i in bb if getattr(self, i, None)]
        dlist += bba # add arguments
        dline = ' '.join(dlist) # to cmd line
        # main args
        cmd = ' '.join([
            '{}'.format(shutil.which('plotHeatmap')),
            '--matrixFile {}'.format(self.matrix),
            '--outFileName {}'.format(self.heatmap_file),
            dline,            
            '1> {}'.format(self.stdout),
            '2> {}'.format(self.stderr),
        ])
        return cmd


    def run(self):
        # save command
        with open(self.heatmap_cmd, 'wt') as w:
            w.write(self.cmd+'\n')
        # run
        if os.path.exists(self.heatmap_file) and not self.overwrite:
            # if re-cal required, remove the old file
            log.info('plotHeatmap() skipped, file exists: {}'.format(self.heatmap_file))
        else:
            log.info('run plotHeatmap: {}'.format(self.matrix))
            os.system(self.cmd)
        # check output
        if not os.path.exists(self.heatmap_file):
            log.error('plotHeatmap() failed, file not found: {}'.format(self.heatmap_file))


def get_args():
    example = ' '.join([
        '$ python bw2heatmap.py', 
        '-m mat.gz -o out_dir',
        '--samplesLabel f1 f2 --regionsLabel g1 g2',
        '--startLabel TSS --endLabel TES -p 8',
    ])
    parser = argparse.ArgumentParser(
        prog='matrix2heatmap.py', description='convert to heatmap', epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser = add_plot_parser(parser)
    # parser.add_argument('-m', dest='matrix', required=True,
    #     help='matrix file, by computeMatrix ')
    # parser.add_argument('-o', dest='out_dir', required=True,
    #     help='directory to save bigWig file')
    # parser.add_argument('-op', '--out-prefix', dest='prefix', default='matrix2heatmap',
    #     help='prefix for output files, default: [matrix2heatmap]')
    # parser.add_argument('--plotType', default='lines',
    #     choices=['lines', 'fill', 'se', 'std'],
    #     help='type of the plot, choices [lines, fill, se, std] default: [lines]')

    # parser.add_argument('-sl', '--samplesLabel', nargs='+', default=None,
    #     help='samples label')
    # parser.add_argument('-rl', '--regionsLabel', nargs='+', default=None,
    #     help='labels for regions in plot, defautl: [None] auto')
    # parser.add_argument('--labelRotation', type=int, default=0, 
    #     help='Rotate the x axis labels in degrees')
    # parser.add_argument('-st', '--startLabel', default='TSS',
    #     help='start label, default: [TSS]')
    # parser.add_argument('-ed', '--endLabel', default='TES',
    #     help='end label, default: [TES]')
    # parser.add_argument('-rf', '--refPointLabel', default='TSS',
    #     help='label on refpoint')

    # parser.add_argument('-cl', '--color-list', dest='colorList', nargs='+', default=None,
    #     help='color list for heatmap, eg: black,yellow,blue default: [None] auto')
    # parser.add_argument('-cm', '--color-map', dest='colorMap', nargs='+', default=None,
    #     help='color map to use for heamtp, eg: Reds Blues')
    # parser.add_argument('-cn', '--color-number', dest='colorNumber', type=int, default=None,
    #     help='control the number of transitions for colorList')
    # parser.add_argument('--alpha', type=float, default=1.0, 
    #     help='the alpha channel, default [1.0]')

    # parser.add_argument('--sortRegions', default='keep',
    #     choices=['descend', 'ascend', 'no', 'keep'],
    #     help='The output should be sorted by the way.')
    # parser.add_argument('--sortUsing', default='mean',
    #     choices=['mean', 'median', 'max', 'min', 'sum', 'region_length'],
    #     help='Which method should be used for sorting, default: [mean]')
    # parser.add_argument('--sortUsingSamples', type=str, default=None,
    #     help='List of sample numbers for sorting, default: [None]')
    # parser.add_argument('--averageTypeSummaryPlot', default='mean',
    #     choices=['mean', 'median', 'max', 'min', 'sum', 'std'],
    #     help='Which method should be used for sorting, default: [mean]')

    # parser.add_argument('--heatmapHeight', type=int, default=12, choices=range(3, 101),
    #     help='Plot height in cm, 3 to 100, default: [12]')
    # parser.add_argument('--heatmapWidth', type=int, default=4, choices=range(1, 101),
    #     help='Plot width in cm, 1 to 100, default: [4]')
    # parser.add_argument('--whatToShow', choices=['plot and heatmap', 'heatmap only', 'heatmap and colorbar'],
    #     default='plot and heatmap', help='plot content, default [plot and heatmap]')

    # parser.add_argument('--yMin', type=float, default=None,
    #     help='Minimum value for the Y-axis')
    # parser.add_argument('--yMax', type=float, default=None,
    #     help='Maximum value for the Y-axis')
    # parser.add_argument('--zMin', type=float, default=None,
    #     help='Minimum value for heatmap intensities')
    # parser.add_argument('--zMax', type=float, default=None,
    #     help='Maximum value for heatmap intensities')
    # parser.add_argument('--perGroup', action='store_true',
    #     help='plot all samples by group')
    # parser.add_argument('-bl', '--blackListFileName', dest='blackListFileName',
    #     default=None, help='blacklist file')
    # parser.add_argument('-p', dest='numberOfProcessors', type=int, default=4,
    #     help='number of processors, default: [4]')
    # parser.add_argument('-O', '--overwrite', dest='overlap', action='store_true',
    #     help='Overwrite output file')
    return parser


def main():
    args = vars(get_args().parse_args())
    Matrix2heatmap(**args).run()

    
if __name__ == '__main__':
    main()

# EOF