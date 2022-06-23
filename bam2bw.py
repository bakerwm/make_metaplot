#!/usr/bin/env python
#-*- encoding:utf-8 -*-
"""
Convert bam to bigwig: bamCoverage

# strand-specific:
fwd: scale(fwd): fwd / total
rev: scale(rev): rev / total

# output
outdir/out.bigWig
outdir/out_fwd.bigWig
outdir/out_rev.bigWig

Example:
$ bamCoverage -b in.bam -o out.bw --binSize 50 \
  --effectiveGenomeSize 100 \
  --scaleFactor 1.0 --normalizeUsing None \
  --blackListFileName b.bed \
  --skipNAs --extendReads --centerReads -p 8
"""

from multiprocessing.sharedctypes import Value
import os
import pathlib
import argparse
import shutil
import pysam
from multiprocessing import Pool
from utils import (
    make_config, update_obj, dump_yaml, file_abspath, file_prefix, symlink_file,
    fix_label, fix_bw, is_valid_file, is_valid_bam, is_valid_bigwig,
    is_valid_bed, log, init_cpu, fix_out_dir, load_yaml
)
from parse_args import get_bam_args, add_bam_parser, add_io_parser
from hiseq.utils.genome import Genome


class Bam2bw(object):
    def __init__(self, **kwargs):
        # c = make_config(**kwargs)
        c = get_bam_args(**kwargs)
        self = update_obj(self, c, force=True)
        self.update_args()


    def update_args(self):
        if isinstance(self.bam_list, str):
            self.bam_list = [self.bam_list] # to list
        if not is_valid_file(self.bam_list, is_valid_bam):
            msg = '\n'.join([
                '{} : {}'.format(is_valid_file(i, is_valid_bam), i) for i in self.bam_list
            ])
            print(msg)
            raise ValueError('bam file illegal')


    def run_single_bam(self, i):
        func = Bam2bw_ss if self.strand_specific else Bam2bw_ns
        args = self.__dict__.copy()
        args.update({'bam': self.bam_list[i], 'prefix': None, })
        return func(**args).run()

        
    def run(self):
        # run in parallel
        if len(self.bam_list) > 1 and self.parallel_jobs > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                out = pool.map(self.run_single_bam, range(len(self.bam_list)))
        else:
            out = [self.run_single_bam(i) for i in range(len(self.bam_list))]
        return out # bw list
        


class Bam2bw_ns(object):
    def __init__(self, **kwargs):
        # c = make_config(**kwargs)
        c = get_bam_args(**kwargs)
        self = update_obj(self, c, force=True)
        self.update_args() # set default to None
        # self.update_labels()
        self.init_files()
        self.cmd = self.get_cmd()
        dump_yaml(self.__dict__, self.config) # save config


    def basic_args(self):
        """
        default arguments [38]
        to-do: outFileSortedRegions, outFileNameMatrix
        """
        # bam, outFileName, outFileFormat
        alist = [
            'bam', 'outFileName',
            'scaleFactor', 'normalizeUsing', 'binSize', 'effectiveGenomeSize',
            'blackListFileName', 'numberOfProcessors',
            'extendReads', 'smoothLength',
            'filterRNAstrand', 'region', 'MNase', 'Offset', 
            'exactScaling', 'ignoreForNormalization',
            'minMappingQuality', 'samFlagInclude', 
            'samFlagExclude', 'minFragmentLength', 'maxFragmentLength'
        ]
        # ['skipNAs', 'centerReads', 'ignoreDuplicates'] # no arguments
        return alist


    def update_args(self):
        alist = self.basic_args()
        d = {i:getattr(self, i, None) for i in alist}
        self = update_obj(self, d, force=True) # update
        if not isinstance(self.prefix, str):
            self.prefix = file_prefix(self.bam, with_path=False)
        if not is_valid_file(self.bam, is_valid_bam):
            raise ValueError('bam file illegal: {}'.format(self.bam))
        # effsize
        es = getattr(self, 'effectiveGenomeSize', None)
        if es is None:
            es = self.get_effsize()
            setattr(self, 'effectiveGenomeSize', es)
        # blacklist
        bl = getattr(self, 'blackListFileName', None)
        if isinstance(self.genome, str):
            if bl is None:
                self.blackListFileName = Genome(self.genome).blacklist()


    def get_effsize(self):
        """
        1. get effective genomesize from default values
        2. return the chromsize from bam file
        """
        effsize = {
            'dm3': 162367812,
            'dm6': 142573017,
            'mm9': 2620345972,
            'mm10': 2652783500,
            'hg19': 2451960000,
            'hg38': 2913022398,
            'GRCh38': 2913022398
        }
        # genome
        genome = getattr(self, 'genome', None)
        s = effsize.get(genome, None)
        # from bam file
        if not isinstance(s, int):
            sam = pysam.AlignmentFile(self.bam)
            s = sum(sam.header.lengths)
        if not isinstance(s, int):
            raise ValueError('could not determin effsize: {}'.format(self.bam))
        return s


    def init_files(self):
        self.out_dir = fix_out_dir(self.out_dir)
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.prefix)
        if not os.path.exists(self.project_dir):
            os.makedirs(self.project_dir)
        self.bam = file_abspath(self.bam)
        prefix = os.path.join(self.project_dir, self.prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.yaml'),
            'bw': prefix+'.bigWig',
            'stdout': prefix+'.bamCoverage.stdout',
            'stderr': prefix+'.bamCoverage.stderr',
            'bw_cmd': prefix+'.bamCoverage.sh',
        }
        self = update_obj(self, args, force=True)
        self.outFileName = self.bw #


    def get_cmd(self):
        """
        construct arguments to command line
        """
        alist = self.basic_args()
        args = {i:getattr(self, i, None) for i in alist}
        dlist = ['--{} {}'.format(k, v) for k,v in args.items() if v is not None]
        bb = ['skipNAs', 'centerReads', 'ignoreDuplicates'] # no arguments
        bba = ['--'+i for i in bb if getattr(self, i, None)]
        dlist += bba # add arguments
        dline = ' '.join(dlist) # to cmd line
        # main args
        cmd = ' '.join([
            '{}'.format(shutil.which('bamCoverage')),
            dline,            
            '1> {}'.format(self.stdout),
            '2> {}'.format(self.stderr),
        ])
        return cmd


    def run(self):
        with open(self.bw_cmd, 'wt') as w:
            w.write(self.cmd+'\n')
        # run
        if os.path.exists(self.bw) and not self.overwrite:
            # if re-cal required, remove the old file
            log.info('bamCoverage() skipped, file exists: {}'.format(self.bw))
        else:
            log.info('run bamCoverage: {}'.format(self.bw))
            os.system(self.cmd)
        # check output
        if not os.path.exists(self.bw):
            log.error('bamCoverage() failed, file not found: {}'.format(self.bw))
        return self.bw


class Bam2bw_ss(object):
    """
    Strand-specific:
    add scale: fwd (fwd/fwd+rev); rev (rev/fwd+rev)
    """
    def __init__(self, **kwargs):
        # c = make_config(**kwargs)
        c = get_bam_args(**kwargs)
        self = update_obj(self, c, force=True)
        self.update_args() # set default to None
        self.init_files()
        dump_yaml(self.__dict__, self.config) # save config


    def update_args(self):
        prefix = getattr(self, 'prefix', None)
        if not isinstance(prefix, str):
            self.prefix = file_prefix(self.bam, with_path=False)
        if not is_valid_file(self.bam, is_valid_bam):
            raise ValueError('bam file illegal: {}'.format(self.bam))


    def init_files(self):
        self.out_dir = fix_out_dir(self.out_dir)
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.prefix)
        self.project_dir = fix_out_dir(self.project_dir)
        prefix = os.path.join(self.project_dir, self.prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.yaml'),
            'bw': prefix+'.bigWig',
            'bw_fwd': prefix+'_fwd.bigWig',
            'bw_rev': prefix+'_rev.bigWig',
        }
        self = update_obj(self, args, force=True)


    def count_bam(self, bam, strand=None):
        """
        for strand-specific RNA-seq (NSR)
        read2 is the sense direction
        -f 16 : forward
        -F 16 : reverse
        """
        c1 = 0 
        c2 = 0
        if strand == '+' or strand is None:
            # (forward) -–filterRNAstrand=forward keeps minus-strand reads, -f 16
            c1 = pysam.view('-c', '-f', '16', '-F', '4', '-@', '8', bam) # fwd
            c1 = int(c1.strip())
        elif strand == '-' or strand is None:
            # (reverse) -–filterRNAstrand=reverse keeps plus-strand reads, -F 16
            c2 = pysam.view('-c', '-F', '16', '-F', '4', '-@', '8', bam) # rev
            c2 = int(c2.strip())
        else:
            log.error('unknown strand: {}'.format(strand))
        return c1 + c2


    def get_bam_count(self, bam, strand=None):
        """
        get BAM count from file, or samtools view -c
        save count to yaml: count.yaml
        """
        bname = file_prefix(bam)
        cf = os.path.join(self.project_dir, bname+'.count.yaml')
        if not os.path.exists(cf):
            cfwd = self.count_bam(bam, strand='+')
            crev = self.count_bam(bam, strand='-')
            d = {
                'name': bname,
                'forward': cfwd,
                'reverse': crev,
                'total': cfwd + crev,
            }
            dump_yaml(d, cf)
        # get from yaml
        d = load_yaml(cf)
        dn = d.get('name', None)
        if strand == '+':
            out = d.get('forward', 0)
        elif strand == '-':
            out = d.get('reverse', 0)
        else:
            out = d.get('total', 0)
        return out


    def get_bam_scale(self, bam):
        """
        the norm scale for forward/reverse strand
        """
        fwd = self.get_bam_count(self.bam, strand='+')
        rev = self.get_bam_count(self.bam, strand='-')
        sf = round(fwd/(fwd+rev), 6)
        sr = round(rev/(fwd+rev), 6)
        return (sf, sr)


    def run(self):
        sf, sr = self.get_bam_scale(self.bam)
        # forward
        args1 = self.__dict__.copy()
        args1.update({
            'out_dir': os.path.join(self.out_dir, self.prefix),
            'prefix': self.prefix+'_fwd',
            'scaleFactor': sf,
            'filterRNAstrand': 'forward',
            # 'normalizeUsing': 'CPM',
        })
        b1 = Bam2bw_ns(**args1)
        b1.run()
        # reverse
        args2 = self.__dict__.copy()
        args2.update({
            'out_dir': os.path.join(self.out_dir, self.prefix),
            'prefix': self.prefix+'_rev',
            'scaleFactor': sf,
            'filterRNAstrand': 'reverse',
            # 'normalizeUsing': 'CPM',
        })
        b2 = Bam2bw_ns(**args2)
        b2.run()
        # link files to upper-level folder
        symlink_file(b1.bw, self.bw_fwd)
        symlink_file(b2.bw, self.bw_rev)
        return (self.bw_fwd, self.bw_rev)


def get_args():
    example = ' '.join([
        '$ bamCoverage -b in.bam -o out.bw --binSize 50',
        '--effectiveGenomeSize 100',
        '--scaleFactor 1.0 --normalizeUsing None',
        '--blackListFileName b.bed',
        '--skkipNAs --extendReads 150 --centerReads -p 8',
    ])
    parser = argparse.ArgumentParser(
        prog='bam2bw.py', description='bam2bw', epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser = add_io_parser(parser)
    parser = add_bam_parser(parser)
    # parser.add_argument('-o', dest='out_dir', required=True,
    #     help='directory to save bigWig file')
    # parser.add_argument('-op', '--out-prefix', dest='prefix', default='bw2matrix',
    #     help='prefix for output files, default: [bw2matrix]')
    # parser.add_argument('-bs', '--binSize', dest='binSize', type=int, default=50,
    #     help='the bin_size, default [50]')
    # parser.add_argument('-g', '--genome', default=None,
    #     help='The reference genome of bam files, default [None]')
    # parser.add_argument('-ss','--strand-specific', dest='strand_specific',
    #     action='store_true', help='Strand-specific, dUTP library')
    # parser.add_argument('-es', '--effsize', dest='effectiveGenomeSize', type=int,
    #     default=None,
    #     help='effective genome size, if not specified, parse from bam header')
    # parser.add_argument('-s', '--scaleFactor', dest='scaleFactor', type=float,
    #     default=1.0,
    #     help='scale factor for the bam, default: [1.0]')
    # parser.add_argument('-n', '--normalizeUsing', default='None',
    #     choices=['RPKM', 'CPM', 'BPM', 'RPGC', 'None'],
    #     help='Use one of the method to normalize reads, default: [None]')
    # parser.add_argument('-bl', '--blackListFileName', dest='blackListFileName',
    #     default=None, help='blacklist file, default: [None]')
    # parser.add_argument('-p', dest='numberOfProcessors', type=int, default=4,
    #     help='number of processors, default: [4]')
    # parser.add_argument('--extendReads', type=int, default=None,
    #     help='extend PE reads to fragment size')
    # parser.add_argument('--centerReads', action='store_true',
    #     help='reads are centered with respect to the fragment length')
    # parser.add_argument('-O', '--overwrite', dest='overlap', action='store_true',
    #     help='Overwrite output file')
    # return parser
    return parser


def main():
    args = vars(get_args().parse_args())
    Bam2bw(**args).run()


if __name__ == '__main__':
    main()














# class Bam2bw(object):
#     def __init__(self, **kwargs):
#         c = make_config(**kwargs)
#         self = update_obj(self, c, force=True)
#         self.init_args()


#     def init_args(self):
#         if isinstance(self.bam, str):
#             self.bam = [self.bam] # to list
#         self.bam = [file_abspath(i) for i in self.bam]
#         if not is_valid_file(self.bam, is_valid_bam):
#             raise ValueError('bam=, illegal, got {}'.format(self.bam[0]))
#         # update out_prefix
#         self.out_prefix = fix_label(self.bam, self.out_prefix)
#         # update parallel jobs
#         self.parallel_jobs = len(self.bam) # n-samples
#         self.numberOfProcessors, self.parallel_jobs = init_cpu(
#             self.numberOfProcessors, self.parallel_jobs
#         )


#     def bam2bw(self, x):
#         args = self.__dict__.copy()
#         i = self.bam.index(x) # get out_prefix
#         args.update({
#             'bam': x,
#             'out_prefix': self.out_prefix[i],
#         })
#         fun = Bam2bw_ss if self.strand_specific else Bam2bw_ns
#         r = fun(**args)
#         return r.run()


#     def run(self):
#         # see help for run parallel in python:
#         # https://www.machinelearningplus.com/python/parallel-processing-python
#         with Pool(processes=self.parallel_jobs) as pool:
#             out = pool.map(self.bam2bw, self.bam) # single args
#             # pool.starmap(self.bam2bw, [[i] for i in self.bam]) # multi args
#         # return [self.bam2bw(i) for i in self.bam]
#         return out


# class Bam2bw_ss(object):
#     """
#     Convert BAM to bigWig, strand-specific, normalizeUsing CPM

#     $ bamCoverage -b in.bam -o out.bw --binSize 50 \
#       --filterRNAstrand reverse/forward

#     default criteria for dUTP library:
#     -–filterRNAstrand=forward keeps minus-strand reads, -f 16
#     -–filterRNAstrand=reverse keeps plus-strand reads, -F 16
#     """
#     def __init__(self, **kwargs):
#         c = make_config(**kwargs)
#         self = update_obj(self, c, force=True)
#         self.init_args()
#         self.init_files()


#     def init_args(self):
#         if isinstance(self.bam, str):
#             tag = is_valid_file(self.bam, is_valid_bam)
#         else:
#             tag = False
#         if not tag:
#             raise ValueError('bam=, expect str, got {}'.format(
#                 type(self.bam).__name__
#             ))
#         self.bam = file_abspath(self.bam)


#     def init_files(self):
#         if not isinstance(self.out_dir, str):
#             self.out_dir = str(pathlib.Path.cwd())
#         self.out_dir = file_abspath(self.out_dir)
#         self.project_dir = os.path.join(self.out_dir, self.out_prefix)
#         prefix = os.path.join(self.project_dir, self.out_prefix)
#         args = {
#             'config': os.path.join(self.project_dir, 'config.yaml'),
#             'bw': prefix+'.bigWig',
#             'bw_fwd': prefix+'_fwd.bigWig',
#             'bw_rev': prefix+'_rev.bigWig',
#             'stdout': prefix+'.bamCoverage.stdout',
#             'stderr': prefix+'.bamCoverage.stderr',
#         }
#         if not os.path.exists(self.project_dir):
#             os.makedirs(self.project_dir)
#         self = update_obj(self, args, force=True)


#     def count_bam_ss(self, x): # strand-specific
#         # split bam by strand: bam_count_dir
#         log.info('counting reads for bam: {}'.format(x))
#         xname = file_prefix(x)
#         c1_file = os.path.join(self.project_dir, xname+'_fwd.count.txt')
#         c2_file = os.path.join(self.project_dir, xname+'_rev.count.txt')
#         # count reads
#         if is_valid_file([c1_file, c2_file], os.path.exists):
#             # (forward)
#             # name,strand,total,count,scale
#             with open(c1_file) as r1:
#                 line = next(r1).strip().split(',') # name,strand,count
#                 s1 = float(line[-1].strip())
#             # (reverse)
#             # name,strand,total,count,scale
#             with open(c2_file) as r2:
#                 line = next(r2).strip().split(',') # name,strand,count
#                 s2 = float(line[-1].strip())
#         else:
#             # (forward) -–filterRNAstrand=forward keeps minus-strand reads, -f 16
#             c1 = pysam.view('-c', '-f', '16', '-F', '4', '-@', '8', x) # fwd
#             c1 = int(c1.strip())
#             # (reverse) -–filterRNAstrand=reverse keeps plus-strand reads, -F 16
#             c2 = pysam.view('-c', '-F', '16', '-F', '4', '-@', '8', x) # fwd
#             c2 = int(c2.strip())
#             # scalefactor
#             s1 = c1/(c1+c2)
#             s2 = c2/(c1+c2)
#             # save to file
#             # name,strand,total,count,scale
#             line1 = '{},{},{},{},{:.6f}'.format(xname, 'forward', c1+c2, c1, s1)
#             line2 = '{},{},{},{},{:.6f}'.format(xname, 'reverse', c1+c2, c2, s2)
#             with open(c1_file, 'wt') as w:
#                 w.write(line1+'\n')
#             with open(c2_file, 'wt') as w:
#                 w.write(line2+'\n')
#         # out: scale
#         return (s1, s2) # fwd, rev


#     def bam2bw_single(self, strand='forward'):
#         # get the scalefactor
#         c1, c2 = self.count_bam_ss(self.bam)
#         s1 = c1/(c1+c2) # forward
#         s2 = c2/(c1+c2) # reverse
#         sf = s2 if strand == 'forward' else s1 #
#         # simplified suffix
#         sx = '_fwd' if strand == 'forward' else '_rev'
#         # run program
#         args = self.__dict__.copy()
#         args.update({
#             'out_dir': os.path.join(self.out_dir, self.out_prefix),
#             'out_prefix': self.out_prefix+sx,
#             'filterRNAstrand': strand,
#             'normalizeUsing': 'CPM',
#             'scaleFactor': sf # scaleFactor
#         })
#         r = Bam2bw_ns(**args)
#         r.run()

#         # copy bw to self.project_dir ?!
#         dest_bw = self.bw_fwd if strand == 'forward' else self.bw_rev
#         symlink_file(r.bw, dest_bw)
#         return r.bw


#     def run(self):
#         return [self.bam2bw_single(i) for i in ['forward', 'reverse']]


# class Bam2bw_ns(object):
#     """ Single BAM file
#     Check: normalizeUsing, scaleFactor, filterRNAstrand (None)
#     --skipNAs
#     --smoothLength
#     --extendReads
#     --centerReads
#     --binSize
#     --effectiveGenomeSize
#     Example:
#     $ bamCoverage -b in.bam -o out.bw --binSize 50 \
#       --effectiveGenomeSize 100 \
#       --scaleFactor 1.0 --normalizeUsing None \
#       --blackListFileName b.bed \
#       --skkipNAs --extendReads --centerReads -p 8
#     """
#     def __init__(self, **kwargs):
#         c = make_config(**kwargs)
#         self = update_obj(self, c, force=True)
#         self.init_args()
#         self.init_files()
#         self.cmd = self.get_cmd()
#         # save config
#         dump_yaml(self.__dict__, self.config)


#     def init_args(self):
#         # check input files
#         if not is_valid_file(self.bam, is_valid_bam):
#             raise ValueError('unable to read BAM file: {}'.format(self.bam))
#         self.bam = file_abspath(self.bam)
#         # effsize
#         self.update_effsize()
#         self.update_blacklist()


#     def get_bam_length(self):
#         s = pysam.AlignmentFile(self.bam)
#         return sum(s.header.lengths)


#     def update_effsize(self):
#         # default effsize
#         s = {
#             'dm3': 162367812,
#             'dm6': 142573017,
#             'mm9': 2620345972,
#             'mm10': 2652783500,
#             'hg19': 2451960000,
#             'hg38': 2913022398,
#             'GRCh38': 2913022398
#         }
#         a = self.effectiveGenomeSize # pre
#         if isinstance(a, int):
#             if a > 0:
#                 out = a
#             else:
#                 out = -1
#         elif self.genome in s:
#             out = s.get(self.genome, -1)
#         else:
#             out = self.get_bam_length()
#         # check
#         if out < 0:
#             raise ValueError('failed, effectiveGenomeSize={}, genome={}'.format(
#                 self.effectiveGenomeSize, self.genome
#             ))
#         self.effectiveGenomeSize = out


#     def update_blacklist(self):
#         if not is_valid_file(self.blacklist, is_valid_bed):
#             self.blacklist = None
#         # self.blacklist = Genome(self.genome).blacklist()
#         # to-do


#     def init_files(self):
#         """
#         set output files
#         """
#         if not isinstance(self.out_dir, str):
#             self.out_dir = str(pathlib.Path.cwd())
#         self.out_dir = file_abspath(self.out_dir)
#         self.project_dir = os.path.join(self.out_dir, self.out_prefix)
#         prefix = os.path.join(self.project_dir, self.out_prefix)
#         args = {
#             'config': os.path.join(self.project_dir, 'config.yaml'),
#             'bw': prefix+'.bigWig',
#             'stdout': prefix+'.bamCoverage.stdout',
#             'stderr': prefix+'.bamCoverage.stderr',
#             'bam2bw_cmd': os.path.join(self.project_dir, prefix+'.bamCoverage.sh')
#         }
#         if not os.path.exists(self.project_dir):
#             os.makedirs(self.project_dir)
#         self = update_obj(self, args, force=True)


#     def get_cmd(self):
#         # basic
#         args_basic = ' '.join([
#             '{}'.format(shutil.which('bamCoverage')),
#             '-b {}'.format(self.bam),
#             '-o {}'.format(self.bw),
#             '-p {}'.format(self.numberOfProcessors),
#             '--skipNAs',
#             '--effectiveGenomeSize {}'.format(self.effectiveGenomeSize),
#         ])
#         if is_valid_file(self.blacklist, is_valid_bed):
#             args_basic += ' -bl {}'.format(self.blacklist)
# #         --extendReads
#         # norm
#         args_norm = ' '.join([
#             '--scaleFactor {}'.format(self.scaleFactor),
#             '--normalizeUsing {}'.format(self.normalizeUsing),
#         ])
#         # strand : filterRNAstrand
#         if hasattr(self, 'filterRNAstrand'):
#             if self.filterRNAstrand in ['forward', 'reverse']:
#                 args_norm += ' --filterRNAstrand {}'.format(self.filterRNAstrand)
#         # binSize:
#         if hasattr(self, 'binSize'):
#             if isinstance(self.binSize, int):
#                 args_norm += ' --binSize {}'.format(self.binSize)
#         # extend
#         if self.extendReads:
#             args_norm += ' --extendReads'
#         # extra
#         args_extra = ' '.join([
#             '--centerReads',
#             '--smoothLength {}'.format(self.binSize*3),
#             '1> {}'.format(self.stdout),
#             '2> {}'.format(self.stderr)
#         ])
#         return ' '.join([args_basic, args_norm, args_extra])


#     def run(self):
#         # save command
#         with open(self.bam2bw_cmd, 'wt') as w:
#             w.write(self.cmd+'\n')
#         # run
#         if os.path.exists(self.bw) and not self.overwrite:
#             # if re-cal required, remove the old file
#             log.info('Bam2bw() skipped, file exists: {}'.format(self.bw))
#         else:
#             log.info('run bamCoverage: {}'.format(self.bw))
#             os.system(self.cmd)
#         # check output
#         if not is_valid_file(self.bw, is_valid_bigwig):
#             log.error('Bam2bw() failed, file failed: {}'.format(self.bw))
#         return self.bw

# EOF