#!/usr/bin/env python
#-*- encoding:utf-8 -*-
"""
Convert bigWig to matrix

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

# for GTF regions
--metagene
--transcriptID transcript
--exonID exon
--transcript_id_designator transcript_id
"""


import os
import re
# import sys
# import pathlib
import argparse
import shutil
from matrix_rbind import Matrix_rbind
from utils import (
    make_config, update_obj, dump_yaml, file_abspath, file_prefix,
    is_valid_bigwig, is_valid_bed, is_valid_gtf, is_valid_file, fix_out_dir,
    load_matrix, log
)
from parse_args import add_io_parser, add_bw_parser, get_bw_args


class Bw2matrix(object):
    def __init__(self, **kwargs):
        # c = make_config(**kwargs)
        c = get_bw_args(**kwargs)
        self = update_obj(self, c, force=True)
        self.update_args()


    def update_args(self):
        # samplesLabel
        if isinstance(self.bw_fwd_list, list) and \
            isinstance(self.bw_rev_list, list):
            if self.samplesLabel is None:
                self.samplesLabel = file_prefix(self.bw_fwd_list)
            self.samplesLabel = [i.replace('_fwd', '') for i in self.samplesLabel]
            self.bw_list = None
            self.strand_specific = True
            self.bw2matrix_func = Bw2matrix_ss
        elif isinstance(self.bw_list, list):
            self.bw_fwd_list, self.bw_rev_list = (None, None)
            if self.samplesLabel is None:
                self.samplesLabel = file_prefix(self.bw_list)
            self.strand_specific = False
            self.bw2matrix_func = Bw2matrix_ns
        else:
            raise ValueError('unknown bigWig files')


    def run(self):
        args = self.__dict__.copy()
        return self.bw2matrix_func(**args).run()


class Bw2matrix_ns(object):
    """
    non-strand-specific
    bw_list:
    bigWig to matrix using: computeMatrix
    """
    def __init__(self, **kwargs):
        # c = make_config(**kwargs)
        c = get_bw_args(**kwargs)
        self = update_obj(self, c, force=True)
        self.update_args() # set default to None
        self.update_labels()
        self.init_files()
        self.cmd = self.get_cmd()
        dump_yaml(self.__dict__, self.config) # save config


    def basic_args(self):
        """
        default arguments [38]
        to-do: outFileSortedRegions, outFileNameMatrix
        no-value-parameters: --metagene, --skipZeros, --missingDataAsZero
        """
        alist = [
            'regionsFileName', 'scoreFileName', 'outFileName',
            'outFileNameMatrix', 'outFileSortedRegions', 'samplesLabel',
            'regionBodyLength', 'unscaled5prime', 'unscaled3prime',
            'referencePoint', 'nanAfterEnd',
            'beforeRegionStartLength', 'afterRegionStartLength', 'binSize',
            'sortRegions', 'sortUsing', 'sortUsingSamples', 'averageTypeBins',
            'minThreshold', 'maxThreshold', 'blackListFileName',
            'scale', 'numberOfProcessors', 'transcriptID', 'exonID',
            'transcript_id_designator'
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
        if self.matrix_type == 'scale-regions':
            self.sub_cmd = 'scale-regions'
            rps = ['referencePoint', 'nanAfterEnd']
            [setattr(self, i, None) for i in rps]
        else:
            self.sub_cmd = 'reference-point'
            srs = [
                'regionBodyLength', 'startLabel', 'endLabel', 'unscaled5prime',
                'unscaled3prime'
            ]
            [setattr(self, i, None) for i in srs]
        if not isinstance(self.prefix, str):
            self.prefix = 'bw2matrix'
        if not is_valid_file(self.bw_list, is_valid_bigwig):
            raise ValueError('bigWig file illegal: {}'.format(self.bw_list))
        # check file extension
        f_ext = [os.path.splitext(i)[1] for i in self.region_list]
        if set(f_ext) == {'.bed'}:
            region_func = is_valid_bed
        elif set(f_ext) == {'.gtf'}:
            region_func = is_valid_gtf
        else:
            raise ValueError('unknown region files:')
        if not is_valid_file(self.region_list, region_func):
            raise ValueError('region files illegal: {}'.format(self.region_list))


    def update_labels(self):
        """
        convert to str, or keep None
        1. samplesLabel
        """
        # 1. samplesLabel
        p = re.compile('^".*"$')
        if isinstance(self.samplesLabel, list):
            k1 = len(self.samplesLabel) == len(self.bw_list) # !!!
            self.samplesLabel = [f'"{i}"' for i in self.samplesLabel if p.match(i) is None] # !!!
            self.samplesLabel = ' '.join(self.samplesLabel)


    def init_files(self):
        self.out_dir = fix_out_dir(self.out_dir)
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = self.out_dir
        # self.project_dir = os.path.join(self.out_dir, self.prefix)
        # if not os.path.exists(self.project_dir):
        #     os.makedirs(self.project_dir)
        # self.bw_list = list(map(lambda i: os.path.abspath(i), self.bw_list))
        self.bw_list = file_abspath(self.bw_list)
        self.region_list = file_abspath(self.region_list)
        prefix = os.path.join(self.project_dir, self.prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.yaml'),
            'matrix': prefix+'.mat.gz',
            'matrix_value': prefix+'.mat.tab',
            'sorted_regions_file': prefix+'.sortedRegions.bed',
            'stdout': prefix+'.computeMatrix.stdout',
            'stderr': prefix+'.computeMatrix.stderr',
            'matrix_cmd': os.path.join(self.project_dir, prefix+'.computeMatrix.sh'),
        }
        self = update_obj(self, args, force=True)


    def get_cmd(self):
        """
        construct arguments to command line
        """
        alist = self.basic_args()
        args = {i:getattr(self, i, None) for i in alist}
        dlist = ['--{} {}'.format(k, v) for k,v in args.items() if v is not None]
        # bb = ['missingDataAsZero', 'skipZeros'] # no arguments
        # bba = ['--'+i for i in bb if getattr(self, i, None)]
        # dlist += bba # add arguments
        for i in ['missingDataAsZero', 'skipZeros', 'metagene']:
            if getattr(self, i, False):
                dlist += ['--'+i]
        dline = ' '.join(dlist) # to cmd line
        # main args
        cmd = ' '.join([
            '{}'.format(shutil.which('computeMatrix')),
            self.sub_cmd,
            '-S {}'.format(' '.join(self.bw_list)),
            '-R {}'.format(' '.join(self.region_list)),
            '-o {}'.format(self.matrix),
            '--outFileNameMatrix {}'.format(self.matrix_value),
            '--outFileSortedRegions {}'.format(self.sorted_regions_file),
            '--numberOfProcessors {}'.format(self.numberOfProcessors),
            dline,
            '1> {}'.format(self.stdout),
            '2> {}'.format(self.stderr),
        ])
        return cmd


    def run(self):
        with open(self.matrix_cmd, 'wt') as w:
            w.write(self.cmd+'\n')
        # run
        if os.path.exists(self.matrix) and not self.overwrite:
            # if re-cal required, remove the old file
            log.info('compureMatrix() skipped, file exists: {}'.format(self.matrix))
        else:
            log.info('run computeMatrix: {}'.format(self.matrix))
            os.system(self.cmd)
        # check output
        if not os.path.exists(self.matrix):
            log.error('computeMatrix() failed, file not found: {}'.format(self.matrix))
        return self.matrix


class Bw2matrix_ss(object):
    """
    strand_specific:
    sense: bw_fwd + region_fwd; bw_rev + region_rev
    anti: bw_fwd + region_rev; bw_rev + region_fwd
    # split region files into strand-specific files
    """
    def __init__(self, **kwargs):
        # c = make_config(**kwargs)
        c = get_bw_args(**kwargs)
        self = update_obj(self, c, force=True)
        self.init_files()
        dump_yaml(self.__dict__, self.config) # save config


    def init_files(self):
        self.out_dir = fix_out_dir(self.out_dir)
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = self.out_dir
        # self.project_dir = os.path.join(self.out_dir, self.prefix)
        self.config = os.path.join(self.project_dir, 'config.yaml')
        prefix = os.path.join(self.project_dir, self.prefix)
        self.matrix_sens = prefix+'_sens.mat.gz'
        self.matrix_anti = prefix+'_anti.mat.gz'


    def split_region_files(self):
        """
        split region files by strand
        """
        # guess file extension
        fext = [os.path.splitext(i)[1] for i in self.region_list]
        if len(set(fext)) > 1:
            raise ValueError('more than 1 type region files found, use BED or GTF only')
        # strand: fwd, +
        # rf_fwd = [os.path.splitext(i)[0] + '_fwd.bed' for i in self.region_list]
        rf_fwd = [os.path.join(self.project_dir, file_prefix(i)+'_fwd'+fext[0]) for i in self.region_list]
        [self.split_region_file(a, b, '+', self.overwrite) for a,b in zip(self.region_list, rf_fwd)]
        # strand: rev, -
        #  rf_rev = [os.path.splitext(i)[0] + '_rev.bed' for i in self.region_list]
        rf_rev = [os.path.join(self.project_dir, file_prefix(i)+'_rev'+fext[0]) for i in self.region_list]
        [self.split_region_file(a, b, '-', self.overwrite) for a,b in zip(self.region_list, rf_rev)]
        # output
        return [rf_fwd, rf_rev] # file lists


    def split_region_file(self, x, out, strand='+', overwrite=False):
        """
        Split region BED/GTF files by strand
        * : original
        + : _fwd
        - : _rev
        """
        xname = file_prefix(x)
        if os.path.exists(out) and overwrite is False:
            print('split_region_file() skipped, file exists: {}'.format(out))
        else:
            ftag = 0
            fext = os.path.splitext(x)[1]
            with open(x) as r, open(out, 'wt') as w:
                for line in r:
                    s = line.strip().split('\t')
                    if len(s) > 5:
                        if fext == '.bed':
                            str_idx = 5
                        elif fext == '.gtf':
                            str_idx = 6
                        else:
                            continue
                        if s[str_idx] == strand or strand == '*':
                            ftag += 1
                            w.write(line)
                    elif strand == '*':
                        ftag += 1
                        w.write(line)
                    # else:
                    #     pass
            # check empty
            if ftag == 0:
                raise ValueError('No region files found, check input file: ', x)
        return out


    def subset_matrix(self, m, o, samples=None, groups=None):
        """
        subset matrix by: samples and groups
        computeMatrixOperations subset -m input.mat.gz -o output.mat.gz --groups "group 1"  --samples "sample 3"
        """
        mh = load_matrix(m, header_only=True)
        sl = mh.get('sample_labels', [])
        rl = mh.get('group_labels', [])
        # update samples
        if isinstance(samples, list) and isinstance(groups, list):
            if all([i in sl for i in samples]) and all([i in rl for i in groups]):
                pass
            else:
                raise ValueError('sampels, groups not found in matrix')
        else:
            raise ValueError('sampels, groups not found in matrix')
        #
        samples = [f'"{i}"' for i in samples]
        groups = [f'"{i}"' for i in groups]
        # print('!A-3', samples)
        cmd = ' '.join([
            '{} subset'.format(shutil.which('computeMatrixOperations')),
            '--samples {}'.format(' '.join(samples)),
            '--groups {}'.format(' '.join(groups)),
            '-m {}'.format(m),
            '-o {}'.format(o),
        ])
        cmd_txt = os.path.join(os.path.dirname(o), file_prefix(o)+'.subset_matrix.sh')
        with open(cmd_txt, 'wt') as w:
                w.write(cmd+'\n')
        # run
        if os.path.exists(o) and not self.overwrite:
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


    def run(self):
        rf_fwd, rf_rev = self.split_region_files()
        # update labels
        sl_fwd = [f'{i}_fwd' for i in self.samplesLabel] # samplesLabel
        sl_rev = [f'{i}_rev' for i in self.samplesLabel] # samplesLabel
        rl_fwd = [os.path.basename(i) for i in rf_fwd] # regionsLabel
        rl_rev = [os.path.basename(i) for i in rf_rev] # regionsLabel
        # run for all
        args = self.__dict__.copy()
        args.update({
            'bw_list': self.bw_fwd_list + self.bw_rev_list,
            'region_list': rf_fwd + rf_rev,
            'samplesLabel': sl_fwd + sl_rev,
            'regionsLabel': rl_fwd + rl_rev,
        })
        # print('!A-1', args['bw_list'])
        b2m = Bw2matrix_ns(**args)
        b2m.run()
        # sense: bw_fwd+region_fwd
        sens1 = os.path.join(self.project_dir, self.prefix+'_sens_1.mat.gz')
        self.subset_matrix(b2m.matrix, sens1, samples=sl_fwd, groups=rl_fwd) # fwd+fwd
        sens2 = os.path.join(self.project_dir, self.prefix+'_sens_2.mat.gz')
        self.subset_matrix(b2m.matrix, sens2, samples=sl_rev, groups=rl_rev) # rev+rev
        # antisense
        anti1 = os.path.join(self.project_dir, self.prefix+'_anti_1.mat.gz')
        anti2 = os.path.join(self.project_dir, self.prefix+'_anti_2.mat.gz')
        self.subset_matrix(b2m.matrix, anti1, samples=sl_fwd, groups=rl_rev) # fwd+rev
        self.subset_matrix(b2m.matrix, anti2, samples=sl_rev, groups=rl_fwd) # rev+fwd
        # 4. merge
        Matrix_rbind(m=[sens1, sens2], o=self.matrix_sens).run()
        Matrix_rbind(m=[anti1, anti2], o=self.matrix_anti).run()
        return [self.matrix_sens, self.matrix_anti]


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
    parser = add_io_parser(parser)
    parser = add_bw_parser(parser)
    # parser.add_argument('-b', dest='bw_list', nargs='+', required=False,
    #     help='bw files')
    # parser.add_argument('-bf', dest='bw_fwd_list', nargs='+', required=False,
    #     help='bw files on forward strand')
    # parser.add_argument('-br', dest='bw_rev_list', nargs='+', required=False,
    #     help='bw files on reverse strand')
    # parser.add_argument('-r', dest='region_list', nargs='+', required=True,
    #     help='region files')
    # parser.add_argument('-o', dest='out_dir', required=True,
    #     help='directory to save bigWig file')
    # parser.add_argument('-op', '--out-prefix', dest='prefix', default='bw2matrix',
    #     help='prefix for output files, default: [bw2matrix]')
    # parser.add_argument('-t', '--matrix-type', dest='matrix_type',
    #     default='scale-regions', choices=['scale-regions', 'reference-point'],
    #     help='choose the matrix type, default: [scale-regions]')
    # parser.add_argument('-sl', '--samplesLabel', nargs='+', default=None,
    #     help='samples label')
    # parser.add_argument('-bs', '--binSize', type=int, default=50,
    #     help='the bin_size, default [50]')
    # parser.add_argument('-u', '--beforeRegionStartLength', type=int, default=500,
    #     help='Distance upstream of TSS, default: [500]')
    # parser.add_argument('-d', '--afterRegionStartLength', type=int, default=500,
    #     help='Distance downstream of TES, default: [500]')
    # parser.add_argument('-m', '--regionBodyLength', type=int, default=1000,
    #     help='Distance for all regions, default: [1000]')
    # parser.add_argument('-st', '--startLabel', default='TSS',
    #     help='start label, default: [TSS]')
    # parser.add_argument('-ed', '--endLabel', default='TES',
    #     help='end label, default: [TES]')
    # parser.add_argument('--sortRegions', default='keep',
    #     choices=['descend', 'ascend', 'no', 'keep'],
    #     help='The output should be sorted by the way.')
    # parser.add_argument('--sortUsing', default='mean',
    #     choices=['mean', 'median', 'max', 'min', 'sum', 'region_length'],
    #     help='Which method should be used for sorting, default: [mean]')
    # parser.add_argument('--sortUsingSamples', type=str, default=None,
    #     help='List of sample numbers for sorting, default: [None]')
    # parser.add_argument('--averageTypeBins',
    #     choices=['mean', 'median', 'min', 'max', 'std', 'sum'],
    #     help='Define the type of method for bin size range, default: [mean]')
    # parser.add_argument('-bl', '--blackListFileName', dest='blackListFileName',
    #     default=None, help='blacklist file')
    # parser.add_argument('-p', dest='numberOfProcessors', type=int, default=4,
    #     help='number of processors, default: [4]')
    # parser.add_argument('-O', '--overwrite', dest='overwrite', action='store_true',
    #     help='Overwrite output file')
    return parser


def main():
    args = vars(get_args().parse_args())
    Bw2matrix(**args).run()


if __name__ == '__main__':
    main()


# EOF