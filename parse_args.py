#!/usr/bin/env python
#-*- encoding:utf-8 -*-
"""
Arguments for metaplot job
"""

import os
import sys
import pathlib
import argparse
import shutil
import logging
from matplotlib import colors


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')


"""
# # bamCoverage
# bam outFileName outFileFormat MNase Offset region exactScaling
# ignoreForNormalization skipNonCoveredRegions ignoreDuplicates
# minMappingQuality  samFlagInclude samFlagExclude minFragmentLength maxFragmentLength

scaleFactor binSize normalizeUsing smoothLength extendReads centerReads
blackListFileName effectiveGenomeSize numberOfProcessors filterRNAstrand 

# # computeMatrix
# regionsFileName scoreFileName outFileName outFileNameMatrix outFileSortedRegions
# smartLabels nanAfterEnd 
# unscaled5prime unscaled3prime sortRegions sortUsing averageTypeBins minThreshold maxThreshold  
# scale  metagene transcriptID exonID transcript_id_designator

binSize regionBodyLength beforeRegionStartLength afterRegionStartLength 
startLabel endLabel referencePoint samplesLabel blackListFileName sortUsingSamples 
numberOfProcessors skipZeros missingDataAsZero

# # plotProfile
# matrixFile outFileName outFileSortedRegions outFileNameData
# hclust silhouette averageType legendLocation plotFileFormat 
kmeans clusterUsingSamples colors numPlotsPerRow plotHeight plotWidth plotType 
startLabel endLabel refPointLabel labelRotation regionsLabel samplesLabel plotTitle 
yAxisLabel yMin yMax perGroup 

# # plotHeatmap
# matrixFile outFileName outFileSortedRegions outFileNameMatrix interpolationMethod 
# dpi silhouette sortRegions sortUsing 
# linesAtTickMarks clusterUsingSamples averageTypeSummaryPlot missingDataColor 
# alpha colorNumber boxAroundHeatmaps plotFileFormat
plotType colorMap colorList whatToShow kmeans sortUsingSamples heatmapHeight heatmapWidth 
zMin zMax yMin yMax xAxisLabel yAxisLabel regionsLabel samplesLabel startLabel endLabel 
refPointLabel labelRotation plotTitle legendLocation perGroup
"""

def get_init_args(**kwargs):
    """
    deeptools 3.5.1
    arguments for bamCoverage, computeMatrix, plotProfile, plotHeatmap
    """
    args = get_bam_args(**kwargs)
    a1 = get_bw_args(**kwargs)
    a2 = get_plot_args(**kwargs)
    args.update(a1)
    args.update(a2)
    return args


def get_bam_args(**kwargs):
    """
    Arguments for bamCoverage of deeptools
    """
    args = {
        'bam': None,
        'binSize': 50,
        'blackListFileName': None,
        'centerReads': True,
        'effectiveGenomeSize': None,
        'extendReads': None,
        'filterRNAstrand': None,
        'genome': None,
        'normalizeUsing': 'None',
        'numberOfProcessors': 12,
        'out_dir': None,
        'overwrite': False,
        'prefix': 'metaplot',
        'scaleFactor': 1.0,
        'smoothLength': 150,
        'strand_specific': False, # for bigWig files, matrix, ...
    }
    args.update(kwargs)
    return args


def get_bw_args(**kwargs):
    """
    Argumetns for computeMatrix of deeptools
    """
    args = {
        'out_dir': None,
        'prefix': None,
        'bw_list': None,
        'bw_fwd_list': None,
        'bw_rev_list': None,
        'samplesLabel': None,
        'region_list': None,
        'startLabel': 'TSS',
        'endLabel': 'TES',
        'beforeRegionStartLength': 2000,
        'regionBodyLength': 5000,
        'afterRegionStartLength': 2000,
        'binSize': 50,
        'referencePoint': 'TSS',
        'matrix_type': 'scale-regions',
        'numberOfProcessors': 12,
        'averageTypeBins': None,
        'blackListFileName': None,
        'missingDataAsZero': True,
        'skipZeros': True,
        'sortUsingSamples': None,
        'strand_specific': False,
        'whatToShow': 'plot, heatmap and colorbar',
    }
    args.update(kwargs)
    return args


def get_plot_args(**kwargs):
    """
    Argumetns for plotProfile, plotHeatmap of deeptools
    """
    args = {     
        'afterRegionStartLength': 2000,
        'averageTypeBins': None,
        'beforeRegionStartLength': 2000,
        'binSize': 50,
        'blackListFileName': None,
        'bw_fwd_list': None,
        'bw_list': None,
        'bw_rev_list': None,
        'clusterUsingSamples': None,
        'colorList': None,
        'colorMap': None,
        'colors': None,
        'endLabel': 'TES',
        'heatmapHeight': 12,
        'heatmapWidth': 4,
        'kmeans': None,
        'labelRotation': None,
        'legendLocation': None,
        'numPlotsPerRow': 2,
        'matrix': None,
        'matrix_type': 'scale-regions',
        'missingDataAsZero': True,
        'numberOfProcessors': 12,
        'out_dir': None,
        'perGroup': True,
        'plotHeight': None,
        'plotTitle': None,
        'plotType': 'lines',
        'plotWidth': None,
        'prefix': None,
        'refPointLabel': 'TSS',
        'referencePoint': 'TSS',
        'regionBodyLength': 5000,
        'regionsLabel': None,
        'region_list': None,
        'samplesLabel': None,
        'skipZeros': True,
        'sortUsingSamples': None,
        'startLabel': 'TSS',
        'strand_specific': False,
        'whatToShow': 'plot, heatmap and colorbar',
        'xAxisLabel': None,
        'yAxisLabel': None,
        'yMin': None,
        'yMax': None,
        'zMin': None,
        'zMax': None,
    }
    args.update(kwargs)
    return args


def add_io_parser(parser):
    if not isinstance(parser, argparse.ArgumentParser):
        log.error('unknown parser')
        return None
    parser.add_argument('-o', dest='out_dir', required=False,
        help='directory to save bigWig file')
    parser.add_argument('-op', '--out-prefix', dest='prefix', default='metaplot',
        help='prefix for output files, default: [metaplot]')
    parser = add_common_parser(parser)
    return parser


def add_common_parser(parser):
    if not isinstance(parser, argparse.ArgumentParser):
        log.error('unknown parser')
        return None    
    parser.add_argument('-bs', '--binSize', dest='binSize', type=int, default=50,
        help='the bin_size, default [50]')
    parser.add_argument('-bl', '--blackListFileName', dest='blackListFileName',
        default=None, help='blacklist file')
    parser.add_argument('-p', dest='numberOfProcessors', type=int, default=4, 
        help='number of processors, default: [4]')
    parser.add_argument('-O', '--overwrite', dest='overwrite', action='store_true',
        help='Overwrite output file')
    parser.add_argument('-sl', '--samplesLabel', nargs='+', default=None,
        help='labels for samples in plot, default: [None] auto')
    parser.add_argument('-st', '--startLabel', default='TSS',
        help='start label, default: [TSS]')
    parser.add_argument('-ed', '--endLabel', default='TES',
        help='end label, default: [TES]')
    return parser


def add_bam_parser(parser):
    if not isinstance(parser, argparse.ArgumentParser):
        log.error('unknown parser')
        return None
    parser.add_argument('-b', dest='bam_list', required=True,
        help='bam files')
    parser.add_argument('-g', '--genome', default=None,
        help='The reference genome of bam files, default [None]')
    parser.add_argument('-ss','--strand-specific', dest='strand_specific',
        action='store_true', help='Strand-specific, dUTP library')
    parser.add_argument('-es', '--effsize', dest='effectiveGenomeSize', type=int,
        default=None,
        help='effective genome size, if not specified, parse from bam header')
    parser.add_argument('-s', '--scaleFactor', dest='scaleFactor', type=float,
        default=1.0,
        help='scale factor for the bam, default: [1.0]')
    parser.add_argument('-n', '--normalizeUsing', default='None',
        choices=['RPKM', 'CPM', 'BPM', 'RPGC', 'None'],
        help='Use one of the method to normalize reads, default: [None]')
    parser.add_argument('--extendReads', type=int, default=None,
        help='extend PE reads to fragment size')
    parser.add_argument('--centerReads', action='store_true',
        help='reads are centered with respect to the fragment length')
    parser.add_argument('-sm', '--smoothLength', type=int, default=None,
        help='smoothlength')
    parser.add_argument('-fs', '--filterRNAstrand', default=None,
        help='filt RNA strand, forward or reverse')
    return parser


def add_bw_parser(parser):
    if not isinstance(parser, argparse.ArgumentParser):
        log.error('unknown parser')
        return None
    parser.add_argument('-bw', dest='bw_list', nargs='+', required=False,
        help='bw files')
    parser.add_argument('-bf', dest='bw_fwd_list', nargs='+', required=False,
        help='bw files on forward strand')
    parser.add_argument('-br', dest='bw_rev_list', nargs='+', required=False,
        help='bw files on reverse strand')
    parser.add_argument('-r', dest='region_list', nargs='+', required=True,
        help='region files')
    parser.add_argument('-t', '--matrix-type', dest='matrix_type',
        default='scale-regions', choices=['scale-regions', 'reference-point'],
        help='choose the matrix type, default: [scale-regions]')
    parser.add_argument('-u', '--beforeRegionStartLength', type=int, default=500,
        help='Distance upstream of TSS, default: [500]')
    parser.add_argument('-d', '--afterRegionStartLength', type=int, default=500,
        help='Distance downstream of TES, default: [500]')
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
    return parser


def add_plot_parser(parser):
    if not isinstance(parser, argparse.ArgumentParser):
        log.error('unknown parser')
        return None
    parser.add_argument('-mx', dest='matrix', required=False,
        help='matrix file, by computeMatrix ')
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
    parser.add_argument('--refPointLabel', default='TSS',
        help='refPointLabel label, default: [TSS]')
    parser.add_argument('--yMin', type=float, default=None,
        help='Minimum value for the Y-axis')
    parser.add_argument('--yMax', type=float, default=None,
        help='Maximum value for the Y-axis')
    parser.add_argument('--perGroup', action='store_true',
        help='plot all samples by group')
    return parser


# # Bam2bw template
# def make_bam2bw_config(**kwargs):
#     # for Bam2bw
#     args_bw = {
#         'bam': None,
#         'out_dir': None,
#         'prefix': 'metaplot',
#         'scaleFactor': 1.0,
#         'normalizeUsing': 'None',
#         'binSize': 100,
#         'numberOfProcessors': 4,
#         'blackListFileName': None,
#         'genome': None,
#         'effectiveGenomeSize': None,
#         'overwrite': False,
#         'strand_specific': False, # for bigWig files, matrix, ...
#         'extendReads': None,
#         'centerReads': False,
#     }
#     args_bw.update(kwargs)
#     return args_bw


# # config template
# def make_config(**kwargs):
#     """
#     Generate config
#     """
#     args_init = {
#         'bam_list': None,
#         'bw_list': None,
#         'bw_fwd_list': None,
#         'bw_rev_list': None,
#         'region_list': None,
#         'out_dir': None,
#         'out_prefix': 'metaplot',
#         'samplesLabel': None, # auto
#         'regionsLabel': None, # auto
#         'colorList': None, # auto
#         'scaleFactor': 1.0,
#         'normalizeUsing': 'None',
#         'binSize': 100,
#         'afterRegionStartLength': 0,
#         'beforeRegionStartLength': 0,
#         'regionBodyLength': 1000,
#         'unscaled5prime': 0,
#         'unscaled3prime': 0,
#         'startLabel': 'TSS',
#         'referencePoint': 'TSS',
#         'refPointLabel': 'TSS',
#         'endLabel': 'TES',
#         'matrix_type': 'scale-regions', # reference-point
#         'numberOfProcessors': 4,
#         'numPlotsPerRow': None,
#         'sortRegions': 'keep',
#         'sortUsing': 'mean',
#         'sortUsingSamples': 1,
#         'averageTypeSummaryPlot': 'mean',
#         'boxAroundHeatmaps': 'yes',
#         'plotTitle': 'Heatmap',
#         'whatToShow': 'plot, heatmap and colorbar',
#         'dpi': 150,
#         'full_version': True,
#         'heatmapHeight': 10,
#         'heatmapWidth': 4,
#         'labelRotation': 0,
#         'yMax': None,
#         'yMin': None,
#         'zMax': None,
#         'zMin': None,
#         'themes': 0,
#         'genome': None,
#         'blacklist': None,
#         'effectiveGenomeSize': None,
#         'averageType': 'mean', # mean, median, min, max, sum and std.
#         'plotType': 'lines', # lines, fill, se, std, overlapped_lines, heatmap
#         'overwrite': False,
#         'strand_specific': False, # for bigWig files, matrix, ...
#         'perGroup': False,
#         'colors': None,
#         'extendReads': None,
#         'matrix': None,
#     }
#     args_init.update(kwargs)
#     return args_init

