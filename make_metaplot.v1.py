#!/usr/bin/env python

"""
Generate metaplot using deeptools toolkit

## Why this script?
In order to generate strand-specific genomebody coverage plots for RNAseq data

    Sense: bigWig(fwd) gene(+); bigWig(rev) gene(-)
Antisense: bigWig(fwd) gene(-); bigWig(fwd) gene(+)

## How-To

### A. DNA seq data (non-stranded)
# eg: ChIP-Seq, ATAC-seq, CUT&Tag, ...
1. BAM -> bigWig
2. bigWig + BED -> matrix
3. matrix -> plots (metaplot, heatmap, ...)

### B. RNA seq data (stranded)
# eg: RNA-Seq, ChrRNAseq, TTseq, ...
1. BAM -> bigWig (fwd, rev)
2. bigWig(fwd) + BED(+) -> matrix-1 (fwd)
3. bigWig(rev) + BED(-) -> matrix-2 (fwd)
4. bigWig(fwd) + BED(-) -> matrix-3 (rev)
5. bigWig(rev) + BED(+) -> matrix-4 (rev)
6. merge matrix (by rows)
   matrix-1 + matrix-2 -> matrix_fwd
   matrix-3 + matrix-4 -> matrix_rev
7. matrix (fwd/rev) -> plots

## Examples

### A. Simple example for DNA
(bam/bw files; region files; out_dir)


### B. Simple example for RNA
(bam files(required); region files; out_dir)


### C. Complicated examples
(normalization, scale, ...)
"""

import os
import sys
import pathlib
import argparse
import shutil
import logging
from matplotlib import colors
from xopen import xopen
import pysam
import pyBigWig
import json
import yaml
import toml

# from hiseq.utils.genome import Genome # !!! ?
# from hiseq.utils.utils import Config, update_obj, log
# from hiseq.utils.file import file_abspath, file_prefix, check_path, file_exists

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')


"""
# Examples
$ computeMatrixOperations subset -m foo.mat.gz -o forward.mat.gz --samples SRR648667.forward SRR648668.forward SRR648669.forward SRR648670.forward
$ computeMatrixOperations subset -m foo.mat.gz -o reverse.mat.gz --samples SRR648667.reverse SRR648668.reverse SRR648669.reverse SRR648670.reverse

# Examples
$ computeMatrixOperations filterStrand -m forward.mat.gz -o forward.subset.mat.gz --strand -
$ computeMatrixOperations filterStrand -m reverse.mat.gz -o reverse.subset.mat.gz --strand +

# Examples
$ computeMatrixOperations rbind -m forward.subset.mat.gz reverse.subset.mat.gz -o merged.mat.gz
$ computeMatrixOperations sort -m merged.mat.gz -o sorted.mat.gz -R genes.gtf
"""

# config template
def make_config(**kwargs):
    """
    Generate config
    """
    args_init = {
        'bam_list': 'None',
        'bw_list': 'None',
        'bw_fwd_list': None,
        'bw_rev_list': None,
        'region_list': 'None',
        'out_dir': 'None',
        'out_prefix': 'metaplot',
        'samplesLabel': 'None', # auto
        'regionsLabel': 'None', # auto
        'colorList': None, # auto
        'scaleFactor': 1.0,
        'normalizeUsing': 'None',
        'binSize': 100,
        'afterRegionStartLength': 0,
        'beforeRegionStartLength': 0,
        'regionBodyLength': 1000,
        'unscaled5prime': 0,
        'unscaled3prime': 0,
        'startLabel': 'TSS',
        'refPointLabel': 'TSS',
        'endLabel': 'TES',
        'matrix_type': 'scale-regions', # reference-point
        'numberOfProcessors': 4,
        'sortRegions': 'keep',
        'sortUsing': 'mean',
        'sortUsingSamples': 1,
        'boxAroundHeatmaps': 'yes',
        'plotTitle': 'Heatmap',
        'whatToShow': 'plot, heatmap and colorbar',
        'dpi': 150,
        'full_version': True,
        'heatmapHeight': 10,
        'heatmapWidth': 4,
        'labelRotation': 0,
        'yMax': 'None',
        'yMin': 'None',
        'zMax': 'None',
        'zMin': 'None',
        'themes': 0,
        'genome': None,
        'blacklist': None,
        'effectiveGenomeSize': None,
        'averageType': 'mean', # mean, median, min, max, sum and std.
        'plotType': 'lines', # lines, fill, se, std, overlapped_lines, heatmap
        'overwrite': False,
        'strand_specific': False, # for bigWig files, matrix, ...
    }
    args_init.update(kwargs)
    return args_init


################################################################################
## 0. help functions
def update_obj(obj, d, force=True, remove=False):
    """Update the object, by dict
    d: dict
    force: bool, update exists attributes
    remove: bool, remove exists attributes
    Update attributes from dict
    force exists attr
    """
    if remove is True:
        for k in obj.__dict__:
            delattr(obj, k)
    # add attributes
    if isinstance(d, dict):
        for k, v in d.items():
            if not hasattr(obj, k) or force:
                setattr(obj, k, v)
    return obj


def load_yaml(x):
    """
    Loding data from YAML file
    x should be file
    """
    d = None
    if not isinstance(x, str):
        return None
    if os.path.exists(x):
        try:
            with open(x, 'r') as r:
                if os.path.getsize(x) > 0:
                    d = yaml.load(r, Loader=yaml.FullLoader)
                    d = dict(sorted(d.items(), key=lambda x:x[0]))
                    # d = collections.OrderedDict(sorted(d.items()))
        except Exception as exc:
            log.error('from_yaml() failed, {}'.format(exc))
        finally:
            return d
    else:
        log.error('from_yaml() failed, file not exists: {}'.format(x))
        

def dump_yaml(d, x):
    """
    Writing data to YAML file
    d dict, data to file
    x str, path to YAML file

    yaml.dump(), does not support OrderedDict
    Solution: OrderedDict -> json -> dict
    """
    # check
    x = os.path.abspath(x)
    if not isinstance(d, dict):
        log.error('to_yaml(d=) failed, dict expect, got {}'.format(
            type(d).__name__))
    elif not isinstance(x, str):
        log.error('to_yaml(d=) failed, str expect, got {}'.format(
            type(x).__name__))
    elif not os.path.exists(os.path.dirname(x)):
        log.error('to_yaml(x=) failed, file not exists: {}'.format(x))
    else:
        try:
            with open(x, 'wt') as w:
                yaml.dump(d, w)
        except:
            log.warning('saving as YOML failed, use TOML instead')
            x_toml = os.path.splitext(x)[0] + '.toml'
            with open(x_toml, 'wt') as w:
                toml.dump(d, w)


def file_abspath(x):
    """Return the absolute path of file
    Parameters
    ----------
    x : str,list
        Path to a file, or list of files
    """
    if x is None or x == 'None': # in case yaml format?!
        out = None
    elif isinstance(x, str):
        out = os.path.abspath(os.path.expanduser(x))
    elif isinstance(x, list):
        out = [file_abspath(i) for i in x]
    else:
        log.warning('x, expect str,list, got {}'.format(type(x).__name__))
        out = x
    return out


def file_prefix(x, with_path=False):
    """Extract the prefix of file,
    compatiabe for None

    Parameters
    ----------
    x : str,list
        Path to a file, or list of files

    remove extensions
    .gz, .fq.gz
    """
    if isinstance(x, str):
        if x.endswith('.gz') or x.endswith('.bz2'):
            x = os.path.splitext(x)[0]
        out = os.path.splitext(x)[0]
        if not with_path:
            out = os.path.basename(out)
    elif isinstance(x, list):
        out = [file_prefix(i, with_path) for i in x]
    elif x is None:
        out = None
    else:
        log.error('unknown x, str,list,None expected, got {}'.format(
            type(x).__name__))
        out = None
    return out


def fix_label(x, labels=None):
    """Make sure the length of samples and labels identical
    support: sample_labels, region_labels, ...

    Parameters
    ----------
    x : list
        Path to the sample, regions

    labels : list
        The labels

    Extract the filename of the samples
    """
    smp_name = [file_prefix(i)[0] for i in x]
    if labels is None:
        labels = smp_name
    if isinstance(labels, str):
        labels = [labels]
    if isinstance(labels, list):
        if not len(x) == len(labels):
            log.error(
                'labels and x not identical,',
                'samples={}, lables={}'.format(len(x), len(labels)))
            labels = smp_name # self.out_filename
        out = labels
    else:
        out = None
    return out


def fix_bw(x):
    """Make sure input files: bigWig exists

    Parameters
    ----------
    x : str or list
        Make sure input bigWig file exists
    """
    if isinstance(x, str):
        # x = file_abspath(x) # absolute path
        try:
            out = is_bw(x)
        except:
            out = False
            log.error('fix_bw() failed, {}'.format(x))
        tag = 'ok' if out else 'failed'
        print('{:>6s} : {}'.format(tag, x))
    elif isinstance(x, list):
        x = [file_abspath(i) for i in x] # absolute path
        out = [self.fix_bw(i) for i in x]
    else:
        log.error('fix_bw(x={}), expect str or list, got {}'.format(x, type(x).__name__))
        out = False
    ## set the absolute the path of the bigwig files
    return out


def is_valid_bam(x):
    """Check if x is valid BAM file

    Parameters
    ----------
    x : str
        Path to the BAM file
    """
    out = False
    if isinstance(x, str):
        if os.path.exists(x):
            try:
                s = pysam.AlignmentFile(x)
                out = True
            except ValueError as err:
                out = False
    return out


def is_valid_bigwig(x):
    """Check if x is bigWig file
    retrieve by keywords
    bw files:
    """
    out = False
    if isinstance(x, str):
        if os.path.exists(x):
            bw = pyBigWig.open(x)
            out = bw.isBigWig()
            bw.close()
    return out


def is_valid_bed(x):
    """Check if x is valid BED file
    file name: .bed, exists
    """
    out = False
    if isinstance(x, str):
        out1 = os.path.exists(x)
        x_ext = os.path.splitext(x)[1]
        out2 = x_ext.lower() in ['.bed', '.bed6', '.bed12', '.narrowpeak']
        out = out1 and out2
    return out


# check_args
# is_valid_bam
# is_valid_bigwig
# is_valid_bed
# is_valid_gtf
# log/print2 # including date
# file_abspath
# file_prefix
# file_exists
# create_dir
#
# strand-specific
#
# deeptools toolkits
# bamCoverage
# computeMatrix
# computeMatrixOperations rbind
# plotProfile
# plotHeatmap



################################################################################
## 1. BAM to BigWig
class Bam2bw(object):
    """ Single BAM file
    Check: normalizeUsing, scaleFactor, filterRNAstrand (None)
    --skipNAs
    --smoothLength
    --extendReads
    --centerReads
    --binSize
    --effectiveGenomeSize
    Example:
    $ bamCoverage -b in.bam -o out.bw --binSize 50 \
      --effectiveGenomeSize 100 \
      --scaleFactor 1.0 --normalizeUsing None \
      --blackListFileName b.bed \
      --skkipNAs --extendReads --centerReads -p 8
    """
    def __init__(self, **kwargs):
        c = make_config(**kwargs)
        self = update_obj(self, c, force=True)
        self.init_args()
        self.init_files()
        self.cmd = self.get_cmd()


    def init_args(self):
        # check input files
        if not is_valid_bam(self.bam):
            raise ValueError('unable to read BAM file: {}'.format(self.bam))
        self.bam = file_abspath(self.bam)
        # effsize
        self.update_effsize()
        self.update_blacklist()


    def get_bam_length(self):
        s = pysam.AlignmentFile(self.bam)
        return sum(s.header.lengths)


    def update_effsize(self):
        a = self.effectiveGenomeSize # pre
        # default effsize
        s = {
            'dm3': 162367812,
            'dm6': 142573017,
            'mm9': 2620345972,
            'mm10': 2652783500,
            'hg19': 2451960000,
            'hg38': 2913022398,
            'GRCh38': 2913022398
        }
        if a in s:
            out = s.get(a, -1)
        elif isinstance(a, int):
            if a > 0:
                out = a
            else:
                out = -1
        elif self.genome in s:
            out = s.get(self.genome, -1)
        else:
            out = self.get_bam_length()
        # check
        if out < 0:
            raise ValueError('unable to get effsize')
        self.effectiveGenomeSize = out


    def update_blacklist(self):
        if not is_valid_bed(self.blacklist):
            self.blacklist = None
        # self.blacklist = Genome(self.genome).blacklist()
        # to-do


    def init_files(self):
        """
        set output files
        """
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.out_prefix)
        prefix = os.path.join(self.project_dir, self.out_prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.bamCoverage.yaml'),
            'bw': prefix+'.bigWig',
            'stdout': prefix+'.bamCoverage.stdout',
            'stderr': prefix+'.bamCoverage.stderr',
            'bam2bw_cmd': os.path.join(self.project_dir, prefix+'.bamCoverage.sh')
        }
        if not os.path.exists(self.project_dir):
            os.makedirs(self.project_dir)
        self = update_obj(self, args, force=True)


    def get_cmd(self):
        # basic
        args_basic = ' '.join([
            '{}'.format(shutil.which('bamCoverage')),
            '-b {}'.format(self.bam),
            '-o {}'.format(self.bw),
            '-p {}'.format(self.numberOfProcessors),
            '--skipNAs',
            '--effectiveGenomeSize {}'.format(self.effectiveGenomeSize),
        ])
        if is_valid_bed(self.blacklist):
            args_basic += ' -bl {}'.format(self.blacklist)
#         --extendReads
        # norm
        args_norm = ' '.join([
            '--scaleFactor {}'.format(self.scaleFactor),
            '--normalizeUsing {}'.format(self.normalizeUsing),
        ])
        # strand : filterRNAstrand
        if hasattr(self, 'filterRNAstrand'):
            if self.filterRNAstrand in ['forward', 'reverse']:
                args_norm += ' --filterRNAstrand {}'.format(self.filterRNAstrand)
        # extra
        args_extra = ' '.join([
            '--centerReads',
            '--smoothLength {}'.format(self.binSize*3),
            '1> {}'.format(self.stdout),
            '2> {}'.format(self.stderr)
        ])
        return ' '.join([args_basic,  args_norm, args_extra])


    def run(self):
        # save command
        with open(self.bam2bw_cmd, 'wt') as w:
            w.write(self.cmd+'\n')
        # run
        if os.path.exists(self.bw) and not self.overwrite:
            # if re-cal required, remove the old file
            log.info('Bam2bw() skipped, file exists: {}'.format(self.bw))
        else:
            log.info('run bamCoverage: {}'.format(self.bw))
            os.system(self.cmd)


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


    def check_files(self):
        c1 = list(map(is_valid_bigwig, self.bw_list))
        c2 = list(map(is_valid_bed, self.region_list))
        msg = '\n'.join([
            '-'*80,
            '>> score files:',
            '\n'.join([
                '{} : {}'.format(i, k) for i,k in zip(c1, self.bw_list)
            ]),
            '>> region files:',
            '\n'.join([
                '{} : {}'.format(i, k) for i,k in zip(c2, self.region_list)
            ]),
            '>> score labels: {}'.format(self.samplesLabel),
            '>> region labels: {}'.format(self.regionsLabel),
            '-'*80,
        ])
        # valid files
        if not all(c1+c2):
            raise Exception('bw_list and region_list, illegal')
        return msg


    def init_files(self):
        """
        set output files
        """
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.out_prefix)
        prefix = os.path.join(self.project_dir, self.out_prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.computeMatrix.yaml'),
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
            '--sortUsingSamples {}'.format(self.sortUsingSamples)
        ])
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


## 3. matrix to plot
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


    def init_args(self):
        # check input files + labels
        d = self.load_header() # matrix exists
        if d is None:
            raise ValueError('unable to read matrix file')
        self.matrix = file_abspath(self.matrix)
        bd = d.get('body', [0])
        rf = d.get('ref point', [None])
        self.is_refpoint = isinstance(rf[0], str) and bd[0] == 0
        # update labels
        d_samples_label = d.get('sample_labels', [None])
        d_regions_label = d.get('group_labels', [None])
        self.is_valid_sl = len(d_samples_label) == len(self.samplesLabel)
        self.is_valid_rl = len(d_regions_label) == len(self.regionsLabel)


    def init_files(self):
        """
        set output files
        """
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.out_prefix)
        prefix = os.path.join(self.project_dir, self.out_prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.plotProfile.yaml'),
            'profile_cmd': os.path.join(self.project_dir, prefix+'.plotProfile.sh'),
            'profile_file': prefix+'.plotProfile.pdf',
            'stdout': prefix+'.plotProfile.stdout',
            'stderr': prefix+'.plotProfile.stderr',
        }
        if not os.path.exists(self.project_dir):
            os.makedirs(self.project_dir)
        self = update_obj(self, args, force=True)


    def load_header(self):
        """
        header in matrix file
        'body':
        'ref point':
        'sample_labels': *
        'group_labels': *
        """
        # read the 1st line from matrix file
        try:
            fh = xopen(self.matrix)
            s = next(fh)
            s = s[1:] # trim the 1st character: "@"
            d = json.loads(s)
        except:
            print('Failed to load matrix file')
            d = None
        return d


    def get_cmd(self):
        # basic
        args_basic = ' '.join([
            '{}'.format(shutil.which('plotProfile')),
            '--matrixFile {}'.format(self.matrix),
            '--outFileName {}'.format(self.profile_file),
            '--dpi {}'.format(self.dpi),
        ])
        if self.is_valid_sl:
            args_basic += ' --samplesLabel {}'.format(' '.join(self.samplesLabel)),
        if self.is_valid_rl:
            args_basic += ' --regionsLabel {}'.format(' '.join(self.regionsLabel)),
        # type
        if self.is_refpoint:
            args_type = '--refPointLabel {}'.format(self.refPointLabel)
        else:
            args_type = ' '.join([
                '--startLabel {}'.format(self.startLabel),
                '--endLabel {}'.format(self.endLabel),
            ])
        # extra
        args_extra = ''
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
        if os.path.exists(self.profile_file) and not self.overwrite:
            # if re-cal required, remove the old file
            log.info('Matrix2profile() skipped, file exists: {}'.format(self.profile_file))
        else:
            log.info('run computeMatrix: {}'.format(self.matrix))
            os.system(self.cmd)


class Matrix2heatmap(object):
    """
    Example:
    $ plotHeatmap \
      -m input.mat -o out.png \
      --whatToShow "plot, heatmap and colorbar" \
      --startLabel TSS --endLabel PAS \ # --refPointLabel TSS
      --heatmapHeight 15 --heatmapWidth 4 \
      --sortRegions descend --sortUsing mean --sortUsingSamples 1 \
      --colorList red,white,blue \
      --samplesLabel WT 15m 30m 120m --regionsLabel PAS-both \
      --plotType lines \ # lines, fill, se, std
      --averageTypeSummaryPlot mean \
      --zMin 0 --zMax 4
    """
    def __init__(self, **kwargs):
        c = make_config(**kwargs)
        self = update_obj(self, c, force=True)
        self.init_args()
        self.init_files()
        self.cmd = self.get_cmd()


    def init_args(self):
        self.matrix = file_abspath(self.matrix)
        # check input files + labels
        d = self.load_header() # matrix exists
        if d is None:
            raise ValueError('unable to read matrix file')
        bd = d.get('body', [0])
        rf = d.get('ref point', [None])
        self.is_refpoint = isinstance(rf[0], str) and bd[0] == 0
        # labels
        d_samples_label = d.get('sample_labels', [None])
        d_regions_label = d.get('group_labels', [None])
        self.is_valid_sl = len(d_samples_label) == len(self.samplesLabel)
        self.is_valid_rl = len(d_regions_label) == len(self.regionsLabel)


    def init_files(self):
        """
        set output files
        """
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.out_prefix)
        prefix = os.path.join(self.project_dir, self.out_prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.plotHeatmap.yaml'),
            'heatmap_cmd': os.path.join(self.project_dir, prefix+'.plotHeatmap.sh'),
            'heatmap_file': prefix+'.plotHeatmap.pdf',
            'stdout': prefix+'.plotHeatmap.stdout',
            'stderr': prefix+'.plotHeatmap.stderr',
        }
        if not os.path.exists(self.project_dir):
            os.makedirs(self.project_dir)
        self = update_obj(self, args, force=True)


    def load_header(self):
        """
        header in matrix file
        'body':
        'ref point':
        'sample_labels': *
        'group_labels': *
        """
        # read the 1st line from matrix file
        try:
            fh = xopen(self.matrix)
            s = next(fh)
            s = s[1:] # trim the 1st character: "@"
            d = json.loads(s)
        except:
            print('Failed to load matrix file')
            d = None
        return d


    def get_cmd(self):
        # basic
        args_basic = ' '.join([
            '{}'.format(shutil.which('plotHeatmap')),
            '--matrixFile {}'.format(self.matrix),
            '--outFileName {}'.format(self.heatmap_file),
            '--dpi {}'.format(self.dpi),
        ])
        if self.is_valid_sl:
            args_basic += ' --samplesLabel {}'.format(' '.join(self.samplesLabel)),
        if self.is_valid_rl:
            args_basic += ' --regionsLabel {}'.format(' '.join(self.regionsLabel)),
        # type
        if self.is_refpoint:
            args_type = '--refPointLabel {}'.format(self.refPointLabel)
        else:
            args_type = ' '.join([
                '--startLabel {}'.format(self.startLabel),
                '--endLabel {}'.format(self.endLabel),
            ])
        # extra
        args_extra = ''
        if isinstance(self.yMin, float):
            args_extra += ' --yMin {}'.format(self.yMin)
        if isinstance(self.yMax, float):
            args_extra += ' --yMax {}'.format(self.yMax)
        if isinstance(self.zMin , float):
            args_extra += ' --zMin {}'.format(self.zMin)
        if isinstance(self.zMax, float):
            args_extra += ' --zMax {}'.format(self.zMax)
        # log
        args_extra += ' 1> {}'.format(self.stdout)
        args_extra += ' 2> {}'.format(self.stderr)
        # cmd
        return ' '.join([args_basic, args_type, args_extra])


    def run(self):
        # save command
        with open(self.heatmap_cmd, 'wt') as w:
            w.write(self.cmd+'\n')
        # run
        if os.path.exists(self.heatmap_file) and not self.overwrite:
            # if re-cal required, remove the old file
            log.info('Matrix2heatmap() skipped, file exists: {}'.format(self.heatmap_file))
        else:
            log.info('run plotHeatmap: {}'.format(self.matrix))
            os.system(self.cmd)


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
        # priority: bw > bam
        if self.is_valid_file(self.bw_list, is_valid_bigwig):
            if isinstance(self.bw_list, str):
                self.bw_list = [self.bw_list]
            self.bw_list = list(map(file_abspath, self.bw_list))
        elif self.is_valid_file(self.bam_list, is_valid_bam):
            if isinstance(self.bam_list, str):
                self.bam_list = [self.bam_list]
            self.bam_list = list(map(file_abspath, self.bam_list))
        else:
            raise ValueError('illegal, bam_list={}, bw_list={}'.format(
                self.bam_list, self.bw_list))
        # check region files (list)
        if self.is_valid_file(self.region_list, is_valid_bed):
            self.region_list = list(map(file_abspath, self.region_list))
        else:
            raise ValueError('region_list, expect list, got {}'.format(
                type(self.region_list).__name__))
        # out files
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        # prefix
        if not isinstance(self.out_prefix, str):
            raise ValueError('out_prefix, expect str, got {}'.format(
                type(self.out_prefix).__name__))
        # labels
        if self.samplesLabel is None or self.samplesLabel == 'None':
            if isinstance(self.bw_list, list):
                self.samplesLabel = file_prefix(self.bw_list) #
            else:
                self.samplesLabel = file_prefix(self.bam_list) #
        if self.regionsLabel is None or self.regionsLabel == 'None':
            self.regionsLabel = file_prefix(self.region_list)


    def is_valid_file(self, x, fun):
        if isinstance(x, str):
            out = fun(x)
        elif isinstance(x, list):
            out = all(list(map(fun, x)))
        else:
            out = False
        return out


    def check_files(self):
        # n_samples
        samples = self.bw_list if self.is_valid_file(self.bw_list, is_valid_bigwig) else self.bam_list
        n_samples = len(samples)
        if not len(samples) == len(self.region_list):
            raise ValueError('samples and lables not identical in length')
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
        """
        set output files
        """
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.out_prefix)
        prefix = os.path.join(self.project_dir, self.out_prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.matrix.yaml'),
            'bam_count_dir': os.path.join(self.project_dir, 'bam_count'),
            'matrix_cmd': os.path.join(self.project_dir, prefix+'.matrix.sh'),
            'matrix': prefix+'.mat.gz',
            'matrix_sens': prefix+'.mat.sens.gz',
            'matrix_anti': prefix+'.mat.anti.gz',
            'stdout': prefix+'.matrix.stdout',
            'stderr': prefix+'.matrix.stderr',
        }
        self = update_obj(self, args, force=True)
        if not os.path.exists(self.bam_count_dir):
            os.makedirs(self.bam_count_dir)

        
    def split_bam_count(self, x):
        ## split bam by strand: bam_count_dir
        xname = file_prefix(x)
        ## -F 16 : forward
        c1_file = os.path.join(self.bam_count_dir, xname+'.count.fwd.txt')
        if os.path.exists(c1_file):
            with open(c1_file) as r:
                line = next(r).strip().split(',') # name,strand,count
                c1 = int(line[2])
        else:
            c1 = pysam.view('-c', '-F', '16', '-F', '4', '-@', '8', x) # fwd
            c1 = int(c1.strip())
            with open(c1_file, 'wt') as w:
                w.write(','.join([xname, 'fwd', str(c1)])+'\n')
        ## -f 16 : reverse 
        c2_file = os.path.join(self.bam_count_dir, xname+'.count.rev.txt')
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
    
    
    def bam2bw(self, x, x_prefix=None): # dna
        args_bam2bw = self.__dict__copy()
        if x_prefix is None:
            x_prefix = file_prefix(x)
        ## (forward)
        args = {
            'bam': x,
            'out_dir': os.path.join(self.out_dir, '01.bam2bw'),
            'out_prefix': x_prefix,
            'normalizeUsing': 'RPGC',
        }
        args_bam2bw.update(args)
        r = Bam2bw(**args_bam2bw)
        r.run()
        return r.bw
    
        
    def split_bam2bw(self, x, x_prefix=None):
        ## a. check norm-scale for fwd+rev # SE reads
        c1, c2 = self.split_bam_count(x)
        args_bam2bw = self.__dict__.copy() 
        if x_prefix is None:
            x_prefix = file_prefix(x)
        ## (forward)
        args_fwd = {
            'bam': x,
            'out_dir': os.path.join(self.out_dir, '01.bam2bw'),
            'out_prefix': x_prefix+'_fwd',
            'filterRNAstrand': 'forward',
            'normalizeUsing': 'CPM', 
            'scaleFactor': c1/(c1+c2) # scaleFactor
        }
        args_bam2bw.update(args_fwd)
        r1f = Bam2bw(**args_bam2bw)
        r1f.run()
        ## (reverse)
        args_rev = {
            'bam': x,
            'out_dir': os.path.join(self.out_dir, '01.bam2bw'),
            'out_prefix': x_prefix+'_rev',
            'filterRNAstrand': 'reverse',
            'normalizeUsing': 'CPM', 
            'scaleFactor': c2/(c1+c2) # scaleFactor
        }
        args_bam2bw.update(args_rev)
        r1r = Bam2bw(**args_bam2bw)
        r1r.run()
        ## output
        return (r1f.bw, r1r.bw)
        
    
    def split_bed(self, x, out_dir, strand='+'):
        """
        Split BED file by strand
        * : original 
        + : _fwd
        - : _rev
        """
        xname = file_prefix(x)
        t = '_rev' if strand == '-' else '_fwd' if strand == '+' else ''
        out = os.path.join(out_dir, xname+t+'.bed')
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
        
        
    def load_header(self, x):
        """
        header in matrix file
        'body':
        'ref point':
        'sample_labels': *
        'group_labels': *
        """
        # read the 1st line from matrix file
        try:
            fh = xopen(x)
            s = next(fh)
            s = s[1:] # trim the 1st character: "@"
            d = json.loads(s)
        except:
            print('Failed to load matrix file')
            d = None
        return d
        
        
    def subset_matrix(self, x, o, samples=None, groups=None):
        """
        Example:
        computeMatrixOperations subset -m input.mat.gz -o output.mat.gz --groups "group 1"  --samples "sample 3"
        """
        d = self.load_header(x)
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
    def run_dna(self):
        """
        For non-stranded BAM files
        """
        # 1. bam to bw
        if not self.is_valid_file(self.bw_list, is_valid_bigwig):
            self.bw_list = [self.bam2bw(i, j) for i,j in zip(self.bam_list, self.samplesLabel)]
        
        # 2. bw to matrix        
        args_bw2matrix = self.__dict__.copy()
        args_tmp = {
            'out_dir': os.path.join(self.out_dir, '02.bw2matrix'),
        }
        args_bw2matrix.update(args_tmp)
        r2 = Bw2matrix(**args_bw2matrix)
        r2.run()
        
        # 3. copy matrix files
        src_dir_rel = os.path.relpath(os.path.dirname(r2.matrix), os.path.dirname(self.matrix))        
        src_rel = os.path.join(src_dir_rel, os.path.basename(r2.matrix))
        os.symlink(src_rel, self.matrix)
        return self.matrix
        
        
    def run_rna(self):
        """
        For stranded BAM files
        """
        # out_dir_fwd, out_dir_rev = [os.path.join(self.out_dir, i) for i in ['fwd', 'rev']] # ?
        # 1. split regions file, by strand
        args_s1 = self.__dict__.copy()
        bed_dir = os.path.join(self.project_dir, 'bed_files')
        if not os.path.exists(bed_dir):
            os.makedirs(bed_dir)
        region_fwd_list = [self.split_bed(i, bed_dir, strand='+') for i in self.region_list]
        region_rev_list = [self.split_bed(i, bed_dir, strand='-') for i in self.region_list]
        
        # 2. bam to bw + matrix
        bw_list = [self.split_bam2bw(i, j) for i,j in zip(self.bam_list, self.samplesLabel)]
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
            'out_dir': os.path.join(self.out_dir, '02.bam2matrix'),
            'out_prefix': self.out_prefix,
            'bw_list': bw_fwd_list + bw_rev_list,
            'region_list': region_fwd_list + region_rev_list,
            'samplesLabel': sl_fwd + sl_rev,
            'regionsLabel': rl_fwd + rl_rev,
        }
        args_bm.update(args_tmp)
        s1 = Bw2matrix(**args_bm)
        s1.run()
        ## subset matrix
        sens1 = os.path.join(args_tmp['out_dir'], self.out_prefix+'.sens_1.mat.gz')
        sens2 = os.path.join(args_tmp['out_dir'], self.out_prefix+'.sens_2.mat.gz')
        anti1 = os.path.join(args_tmp['out_dir'], self.out_prefix+'.anti_1.mat.gz')
        anti2 = os.path.join(args_tmp['out_dir'], self.out_prefix+'.anti_2.mat.gz')
        self.subset_matrix(s1.matrix, sens1, samples=sl_fwd, groups=rl_fwd) # fwd+fwd
        self.subset_matrix(s1.matrix, sens2, samples=sl_rev, groups=rl_rev) # rev+rev
        self.subset_matrix(s1.matrix, anti1, samples=sl_fwd, groups=rl_rev) # fwd+rev
        self.subset_matrix(s1.matrix, anti2, samples=sl_rev, groups=rl_fwd) # rev+fwd
        ## merge
        self.rbind_matrix(x=[sens1, sens2], o=self.matrix_sens)
        self.rbind_matrix(x=[anti1, anti2], o=self.matrix_anti)
        return (self.matrix_sens, self.matrix_anti)
        
        
    def run(self):
        print(self.check_files())
        if self.strand_specific:
            self.run_rna()
        else:
            self.run_dna()
            

class Bam2profile(object):
    """
    Convert BAM to profile: strand-specific
    """    
    def __init__(self, **kwargs):
        c = make_config(**kwargs)
        self = update_obj(self, c, force=True)
        self.init_args()
        self.init_files()
        self.check_files()
        
    
    def init_args(self):
        pass
    
    
    def check_files(self):
        # config file
        if not isinstance(self.config, str):
            raise ValueError('illegal config={}'.format(self.config))
        config_dir = os.path.dirname(self.config)
        if not os.path.exists(config_dir):
            os.makedirs(config_dir)
    

    def init_files(self):
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.out_prefix)
        prefix = os.path.join(self.project_dir, self.out_prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.yaml'),
            'profile_file': prefix+'.profile.pdf',
            'profile_file_sens': prefix+'.profile_sens.pdf',
            'profile_file_anti': prefix+'.profile_anti.pdf',
            'stdout': prefix+'.profile.stdout',
            'stderr': prefix+'.profile.stderr',
        }
        self = update_obj(self, args, force=True)
        if not os.path.exists(self.project_dir):
            os.makedirs(self.project_dir)        
        

    def run(self):
        # 1. bam to matrix
        args1 = self.__dict__.copy()
        args_tmp = {
            'out_dir': os.path.join(self.out_dir, '01.bam2matrix'),
        }
        args1.update(args_tmp)
        r1 = Bam2matrix(**args1)
        r1.run()
        
        # 2. matrix to profile
        if self.strand_specific:
            # sens
            args2 = self.__dict__.copy()
            args2.update({
                'matrix': r1.matrix_sens,
                'out_prefix': self.out_prefix+'_sens',
                'out_dir': os.path.join(self.out_dir, '02.matrix2profile'),
            })
            r2 = Matrix2profile(**args2)
            r2.run()
            
            # anti
            args3 = self.__dict__.copy()
            args3.update({
                'matrix': r1.matrix_anti,
                'out_prefix': self.out_prefix+'_anti',
                'out_dir': os.path.join(self.out_dir, '02.matrix2profile'),
            })
            r3 = Matrix2profile(**args3)
            r3.run()
            # output
            p = [r2.profile_file, r3.profile_file]
            
            # 3. copy files
            if not os.path.exists(self.profile_file_sens) or self.overwrite:
                shutil.copy(r2.profile_file_sens, self.profile_file_sens)
            if not os.path.exists(self.profile_file_anti) or self.overwrite:
                shutil.copy(r2.profile_file_anti, self.profile_file_anti)
        else:
            args2.update(args_tmp)
            r2 = Matrix2profile(**args2)
            r2.run()
            p = r2.profile_file
            
            # 3. copy files
            if os.path.exists(self.profile_file) and not self.overwrite:
                log.info('bam2profile() skipped, file exists: {}'.format(
                    self.profile_file))
            else:
                log.info('run bam2profile(): {}'.format(self.profile_file))
                shutil.copy(r2.profile_file, self.profile_file)


class Bw2profile(object):
    """
    Convert BigWig to profile
    require: bw_fwd_list, bw_rev_list, bw_list
    """    
    def __init__(self, **kwargs):
        c = make_config(**kwargs)
        self = update_obj(self, c, force=True)
        self.init_args()
        self.init_files()
        self.check_files()


    def init_args(self):
        pass


    def check_files(self):
        if not isinstance(self.config, str):
            raise ValueError('illegal config={}'.format(self.config))
        config_dir = os.path.dirname(self.config)
        if not os.path.exists(config_dir):
            os.makedirs(config_dir)
            
            
    def is_valid_file(self, x, fun):
        if isinstance(x, str):
            out = fun(x)
        elif isinstance(x, list):
            out = all(list(map(fun, x)))
        else:
            out = False
        return out

    
    def check_bw_files(self):
        """
        Strand-specific: True/False
        """
        if all([self.is_valid_file(i, is_valid_bigwig) for i in self.bw_list]):
            self.strand_specific = False
        elif all([self.is_valid_file(i, is_valid_bigwig) for i in self.bw_fwd_list+self.bw_rev_list]):
            self.strand_specific = True
            if len(self.bw_fwd_list) == len(self.bw_rev_list):
                for f,r in zip(self.bw_fwd_list, self.bw_rev_list):
                    f1 = file_prefix(f).strip('_sens')
                    r1 = file_prefix(r).strip('_anti')
                    if not f1 == r1:
                        raise ValueError('bw not match: {} {}'.format(f, r))
            else:
                raise ValueError('bw not match: {} {}'.format(
                    len(self.bw_fwd_list), len(self.bw_rev_list)
                ))
        else:
            raise ValueError('bw_list, bw_fwd_list, bw_rev_list required:')
            
            
    def init_files(self):
        if not isinstance(self.out_dir, str):
            self.out_dir = str(pathlib.Path.cwd())
        self.out_dir = file_abspath(self.out_dir)
        self.project_dir = os.path.join(self.out_dir, self.out_prefix)
        prefix = os.path.join(self.project_dir, self.out_prefix)
        args = {
            'config': os.path.join(self.project_dir, 'config.yaml'),
            'profile_file': prefix+'.profile.pdf',
            'profile_file_sens': prefix+'.profile_sens.pdf',
            'profile_file_anti': prefix+'.profile_anti.pdf',
            'stdout': prefix+'.profile.stdout',
            'stderr': prefix+'.profile.stderr',
        }
        self = update_obj(self, args, force=True)
        if not os.path.exists(self.project_dir):
            os.makedirs(self.project_dir)

            
    def split_bed(self, x, out_dir, strand='+'):
        """
        Split BED file by strand
        * : original 
        + : _fwd
        - : _rev
        """
        xname = file_prefix(x)
        t = '_rev' if strand == '-' else '_fwd' if strand == '+' else ''
        out = os.path.join(out_dir, xname+t+'.bed')
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
            
    
    def bw2matrix_ss(self):
        """
        Strand-specific: True
        """
        # 1. split regions file, by strand
        args_s1 = self.__dict__.copy()
        bed_dir = os.path.join(self.project_dir, '01.bed_files')
        if not os.path.exists(bed_dir):
            os.makedirs(bed_dir)
        region_fwd_list = [self.split_bed(i, bed_dir, strand='+') for i in self.region_list]
        region_rev_list = [self.split_bed(i, bed_dir, strand='-') for i in self.region_list]
        # labels
        sl_fwd = [i+'_fwd' for i in self.samplesLabel]
        sl_rev = [i+'_rev' for i in self.samplesLabel]
        rl_fwd = [os.path.basename(i) for i in region_fwd_list] # 
        rl_rev = [os.path.basename(i) for i in region_rev_list] #

        # 2. bw to matrix
        args2 = self.__dict__.copy()
        args2.update({
            'out_dir': os.path.join(self.out_dir, '02.bw2matrix'),
            'out_prefix': self.out_prefix,
            'bw_list': self.bw_fwd_list + self.bw_rev_list,
            'region_list': region_fwd_list + region_rev_list,
            'samplesLabel': sl_fwd + sl_rev,
            'regionsLabel': rl_fwd + rl_rev,
        })
        r2 = Bw2matrix(**args2)
        r2.run()
        
        # 3. subset matrix
        sens1 = os.path.join(args2['out_dir'], self.out_prefix+'.sens_1.mat.gz')
        sens2 = os.path.join(args2['out_dir'], self.out_prefix+'.sens_2.mat.gz')
        anti1 = os.path.join(args2['out_dir'], self.out_prefix+'.anti_1.mat.gz')
        anti2 = os.path.join(args2['out_dir'], self.out_prefix+'.anti_2.mat.gz')
        self.subset_matrix(r2.matrix, sens1, samples=sl_fwd, groups=rl_fwd) # fwd+fwd
        self.subset_matrix(r2.matrix, sens2, samples=sl_rev, groups=rl_rev) # rev+rev
        self.subset_matrix(r2.matrix, anti1, samples=sl_fwd, groups=rl_rev) # fwd+rev
        self.subset_matrix(r2.matrix, anti2, samples=sl_rev, groups=rl_fwd) # rev+fwd
        
        # 4. merge
        matrix_sens = os.path.join(args2['out_dir'], self.out_prefix+'.sens.mat.gz')
        matrix_anti = os.path.join(args2['out_dir'], self.out_prefix+'.anti.mat.gz')
        self.rbind_matrix(x=[sens1, sens2], o=matrix_sens)
        self.rbind_matrix(x=[anti1, anti2], o=matrix_anti)
        return (matrix_sens, matrix_anti)
    

    def bw2matrix_ns(self):
        """
        Strand-specific: False
        """
        # 1. bw to matrix
        args2 = self.__dict__.copy()
        args2.update({
            'out_dir': os.path.join(self.out_dir, '02.bw2matrix'),
            'out_prefix': self.out_prefix,
            'bw_list': self.bw_list,
            'region_list': self.region_list,
            'samplesLabel': self.samplesLabel,
            'regionsLabel': self.regionsLabel,
        })
        r2 = Bw2matrix(**args2)
        r2.run()
        return r2.matrix

    
    def run(self):
        self.check_bw_files()
        if self.strand_specific:
            # 1. bw to matrix
            matrix_sens, matrix_anti = self.bw2matrix_ss()
            # 2. matrix to profile: sens
            args2 = self.__dict__.copy()
            args2.update({
                'matrix': matrix_sens,
                'out_dir': os.path.join(self.out_dir, '03.matrix2profile'),
            })
            r2 = Matrix2profile(**args2)
            r2.run()
            # 3. copy files: sens
            if not os.path.exists(self.profile_file_sens) or self.overwrite:
                shutil.copy(r2.profile_file, self.profile_file_sens)
            # 4 matrix to profile: anti
            args3 = self.__dict__.copy()
            args3.update({
                'matrix': matrix_anti,
                'out_dir': os.path.join(self.out_dir, '03.matrix2profile'),
            })
            r3 = Matrix2profile(**args3)
            r3.run()
            # 5 copy files: anti
            if not os.path.exists(self.profile_file_anti) or self.overwrite:
                shutil.copy(r2.profile_file, self.profile_file_anti)
        else:
            # 1. bw to matrix
            matrix = self.bw2matrix_ns()
            # 2. matrix to profile
            args2 = self.__dict__.copy()
            args2.update({
                'matrix': matrix,
                'out_dir': os.path.join(self.out_dir, '03.matrix2profile'),
            })
            r2 = Matrix2profile(**args2)
            r2.run()
            # 3. copy files
            if not os.path.exists(self.profile_file) or self.overwrite:
                shutil.copy(r2.profile_file, self.profile_file)

    

def get_args():
    """
    Parsing arguments for plotHeatmaps, deeptools
    """
    example = '\n'.join([
        'Example:',
        '# 1. Generate template config file',
        '$ python make_metaplot.py -c a.yaml -t',
        '# 2: Run program',
        '# modify the config file `a.yaml` according to your data',
        '$ python make_metaplot.py -c a.yaml',
    ])
    parser = argparse.ArgumentParser(prog='make_metaplot',
                                     description='plotProfile',
                                     epilog=example,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--config', default=None, required=True,
        help='configs in .yaml file')
    parser.add_argument('-t', '--get-template', dest='get_template', action='store_true',
        help='Generate the template arguments')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite the plot files')
    return parser




def main():
    args = vars(get_args().parse_args())
#     Bam2profile(**args)
    if args['get_template']:
        args_t = make_config() #
        dump_yaml(args_t, args['config'])
        # log
        msg = '\n'.join([
            '-'*80,
            '# 1. Generating the template config file ...',
            '$ python make_profile.py -t -c {}'.format(args['config']),
            '# 2. Modify the values in YAML' ,
            '# Attentation to the following fields:',
            '  - bam_list: '
            '  - bw_list: ',
            '  - region_list: ',
            '  - samplesLabel: a b c',
            '  - regionsLabel: e f g',
            '  - colorList: red,white,blue',
            '  - yMin, yMax',
            '# 3. Run the command again',
            '$ python make_profile.py -c {}'.format(args['config']),
            '-'*80,
        ])
        print(msg)
    else:
        args_t = make_config()
        args_c = load_yaml(args['config'])
        if args_c is None:
            print('unable to read config: {}'.format(args['config']))
        else:
            args_t.update(args_c)
            r = Bam2profile(**args_t)
            r.run()
    
    
if __name__ == '__main__':
    main()

