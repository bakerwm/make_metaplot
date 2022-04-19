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
        'bam_list': None,
        'bw_list': None,
        'bw_fwd_list': None,
        'bw_rev_list': None,
        'region_list': None,
        'out_dir': None,
        'out_prefix': None,
        'samplesLabel': None, # auto
        'regionsLabel': None, # auto
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
        'perGroup': False,
        'colors': None,
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


def init_cpu(threads=1, parallel_jobs=1):
    """
    Check threads, parallel_jobs
    """
    n_cpu = os.cpu_count() # alternative: multiprocessing.cpu_count()
    max_jobs = int(0.8 * n_cpu) # max 80% cups
    if parallel_jobs > max_jobs:
        parallel_jobs = max_jobs
    if threads > max_jobs:
        threads = max_jobs
    # check parallel_jobs (max: 1/2 of n_cpus)
    # priority: parallel_jobs > threads
    if parallel_jobs*threads > max_jobs:
        # threads [1, max_jobs]
        # parallel_jobs [1, max_jobs]
        threads = int(max_jobs / parallel_jobs)
#     if parallel_jobs > max_jobs: 
#         log.warning('Too large, change parallel_jobs from {} to {}'.format(
#             parallel_jobs, max_jobs))
#         parallel_jobs = max_jobs
#     ## check threads
#     max_threads = int(0.8 * n_cpu / parallel_jobs)
#     if threads * parallel_jobs > 0.8 * n_cpu:
#         log.warning('Too large, change threads from {} to {}'.format(
#             threads, max_threads))
#         threads = max_threads
    return (threads, parallel_jobs)


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


def symlink_file(src, dest, absolute_path=False, force=False):
    """Create symlink, 
    
    Parameters
    ----------
    src : str
        The source file
        
    dest : str
        The destinate file, dir exists
        
    absolute_path : bool
        Use abs_path instead

    force : bool
        Copy files, overwrite dest file
    """
    if not isinstance(src, str):
        log.error('src, expect str, got {}'.format(type(src).__name__))
    elif not isinstance(dest, str):
        log.error('dest, expect str, got {}'.format(type(src).__name__))
    elif os.path.isfile(src):
        src = os.path.abspath(os.path.expanduser(os.path.expandvars(src)))
        if os.path.isdir(dest):
            dest_file = os.path.join(dest, os.path.basename(src))
        else:
            dest_file = dest
        dest_file = os.path.abspath(os.path.expanduser(os.path.expandvars(dest_file)))
        # the relative path of src
        src_dir = os.path.dirname(src)
        src_name = os.path.basename(src)
        dest_dir = os.path.dirname(dest_file)
        src_dir_rel = os.path.relpath(src_dir, dest_dir)
        src_rel = os.path.join(src_dir_rel, src_name)
        src_file = src if absolute_path else src_rel
        # do-the-thing
        if os.path.exists(dest_file) and not force:
            log.info('symlink_file() skipped, dest exists: {}'.format(dest_file))
        else:
            try:
                os.symlink(src_file, dest_file)
            except:
                log.error('symlink_file() failed, {}'.format(dest_file))
    elif os.path.islink(src):
        pass
    else:
        log.warning('symlink_file() failed, src not vaild: {}'.format(src))


def fix_label(x, labels=None):
    """
    Make sure the length of samples and labels identical
    support: sample_labels, region_labels, ...

    Parameters
    ----------
    x : list
        Path to the sample, regions

    labels : list
        The labels

    Extract the filename of the samples
    """
    if isinstance(x, str):
        x = [x]
    if isinstance(labels, str):
        labels = [labels]
    smp_name = [file_prefix(i) for i in x]
    if isinstance(labels, list):
        out = labels if len(x) == len(labels) else smp_name
    else:
        out = smp_name
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


def is_valid_file(x, fun):
    """
    Check if x is valid file: BAM|bigWig|BED

    Parameters
    ----------
    x : str
        Path to the file
    
    fun : function
        The function, is_valid_bam, is_valid_bigwig, is_valid_bed
    """
    if isinstance(x, str):
        out = fun(x)
    elif isinstance(x, list):
        out = all(list(map(fun, x)))
    else:
        out = False
    return out


def is_valid_bam(x):
    """
    Check if x is valid BAM file

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
    """
    Check if x is bigWig file
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
    """
    Check if x is valid BED file
    file name: .bed, exists
    """
    out = False
    if isinstance(x, str):
        out1 = os.path.exists(x)
        x_ext = os.path.splitext(x)[1]
        out2 = x_ext.lower() in ['.bed', '.bed6', '.bed12', '.narrowpeak']
        out = out1 and out2
    return out


def load_matrix_header(x):
    """
    header in matrix file (gzipped)
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
        print('Failed to load matrix file: {}'.format(x))
        d = None
    return d
