#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Extract all SNP from specific region: chr1:100-200 (no more than 10 kb)
Input: 
  - Alignment: in BAM format, sorted 
  - Reference: in FASTA format
  - region: string "chr1:100-200" or in BED format
To-Do: 
  - DNA or RNA mode
commands in details: 
  - deeptools mpileup 
Alignment: BWA with default parameters
bwa mem -t ${CPU} ${index} ${fq}
DNA: strandless 
RNA: strand +/- [forward: -f 16; reverse: -F 16]
pipeline
1. samtools view: separate forward strand (-f 16) and reverse strand (-F 16)
2. samtools mpileup: filt and extract nucleotide content at each position
3. sequenza-utils pileup2acgt: extract ACGT content and frequency
4. functions: further filt editsites
output:
1       chrom
2       chromStart
3       chromEnd
4       name of edits, [chrom]_[chromEnd]_[read_depth]_[pct_of_edit]%
5       REF
6       ALT
7       read_depth
8       A count
9       C count
10      G count
11      T count
12      N count       
"""

__author__ = "Ming Wang <wangm08@hotmail.com>"
__copyright__ = "2019 by Ming Wang <wangm08@hotmail.com>"
__license__ = "MIT"
__email__ = "wangm08@hotmail.com"
__version__ = "0.1"


import os
import sys
import argparse
import tempfile
import logging
import shlex
import subprocess
import logging
import pathlib
import pysam


logging.basicConfig(format='[%(asctime)s] %(message)s', 
                    datefmt='%Y-%m-%d %H:%M:%S', 
                    level=logging.DEBUG)


def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog = 'parse_edits', description='Parsing editing sites',
        epilog='python parse_edits.py -i test.bam -f ref.fa -o edits')
    parser.add_argument('-i', '--bam-list', dest='bam_list', nargs='+', required=True,
        help='BAM files')
    parser.add_argument('-f', '--ref', required=True, 
        help='Reference sequence, FASTA format')
    parser.add_argument('-o', '--out-dir', dest='out_dir', required=True, 
        help='file to save the results')
    parser.add_argument('-n', '--smp-name', nargs='+', dest='smp_name', default=None,
        help='Name of the sample')
    parser.add_argument('-d', '--max-depth', dest='max_depth', default=250,
        type=int, 
        help='maximum read depth at editing position, default: [250]')
    parser.add_argument('-c', '--cutoff', dest='cutoff', default=0.1, 
        type=float, 
        help='minimum editing percentage [0-1], default: [0.1]')
    parser.add_argument('-r', '--region', default=None,
        help='Comma separated list of regions in which pileup is generated')
    parser.add_argument('-R', '--region-list', dest='region_list', default=None,
        help='FILE of region list')
    # parser.add_argument('-t', default='DNA', choices=['DNA', 'RNA'],
    #     help='call SNP for RNA or DNA, default: DNA')
    parser.add_argument('--overwrite', action='store_true',
        help='if specified, overwrite existing file')
    # args = parser.parse_args()
    return parser


# 1. extract all variants
# input: bcftools mpileup -Ob -a FORMAT/DP,FORMAT/AD -f ref.fa test.bam > test.mpileup.bcf
# output: 
# by pysam

"""
DP: Read Depth, refers to the overall read depth from all target samples supporting the genotype call.
DP4: field have four subfields for sequence reads covering the variant. These subfields refer to 
    reference allele covered by forward read, reference allele covered by reverse read, alternate 
    allele covered by forward read, and alternate allele covered by reverse read. DP4 may not sum 
    to DP as it excludes low-quality bases.
AD: Allelic Depths for the ref and alt alleles in the order listed
    AD refers to the allele depth. AD reports the informative reads supporting each allele. 
    AD may not always sum to DP.
"""

def run_mpileup_single(bam, ref, **kwargs):
    """
    Parameters:
    -----------
    bam : str
        path to a bam file, index-ed
    ref : str
        path to the reference file in fasta format
    optional:
    ---------
    out_dir : str
        save the output files, default [pwd] 
    smp_name : str
        prefix of the output files, default [basename(bam)]
    max_depth : int 
        maximum depth allowed, default [250] 
    region : str
        regions in which pileup is generated, comma separated 
    region_list : str
        FILE with regions 
    """
    args = {
        'out_dir': None,
        'smp_name': None,
        'max_depth': 250,
        'threads': 8,
        'region': None,
        'region_file': None,
        'overwrite': False,
    }
    args.update(kwargs)
    # prepare dir    
    if not isinstance(args['out_dir'], str):
        args['out_dir'] = pathlib.Path.cwd()
    args['out_dir'] = os.path.abspath(args['out_dir'])
    if not os.path.exists(args['out_dir']):
        os.makedirs(args['out_dir'])
    # prepare file
    bam = os.path.abspath(bam)
    ref = os.path.abspath(ref)
    if not isinstance(args['smp_name'], str):
        args['smp_name'] = os.path.basename(os.path.splitext(bam)[0])
    bcf_out = os.path.join(args['out_dir'], args['smp_name']+'.mpileup.bcf')
    # basic
    cmd = ' '.join([
        'bcftools mpileup -a FORMAT/DP,FORMAT/AD -Ob',
        f'-d {args["max_depth"]}',
        f'-f {ref}',
        f'-o {bcf_out}',
        f'--threads {args["threads"]}',
        bam,
    ])
    # extra: region
    if isinstance(args['region_list'], str):
        cmd += f' -R {args["region_list"]}'
    elif isinstance(args['region'], str):
        cmd += f' -r {args["region"]}'
    else:
        pass
    # extra: index 
    cmd += f' && bcftools index {bcf_out}'
    # save cmd
    cmd_sh = os.path.join(args['out_dir'], args['smp_name']+'.run_mpileup.sh')
    with open(cmd_sh, 'wt') as w:
        w.write(cmd+'\n')
    # run
    try:
        if os.path.exists(bcf_out) and not args['overwrite']:
            print('parse_edits ... skipped, file exists')
        else:
            os.system(cmd)
    except:
        print('Failed')
    return bcf_out


def run_mpileup(bam_list, ref, **kwargs):
    # fix input
    if not isinstance(bam_list, list):
        print(f'bam_list, require list, got {type(bam_list).__name__}')
        return None
    # fix smp_name
    if isinstance(kwargs['smp_name'], list):
        if len(kwargs['smp_name']) == len(bam_list):
            pass
        
    print(f'Processing bam files: ')
    n_bam = len(bam_list)
    i = 0
    for bam in bam_list:
        i += 1
        print(f'- [{i:>2d}/{n_bam}] {os.path.basename(bam)}')
        bcf = run_mpileup_single(bam, ref, **kwargs)
        parse_mpileup(bcf)


# only for the first sample
def read_mpileup(x):
    """
    Read the bcftools mpileup output, required AD and DP format
    1. run bcftools (version 1.16): 
    $ bcftools mpileup -a FORMAT/DP,FORMAT/AD -f ref.fa -Ob test.bam 
    2. extract base content at each position (specified, no more than 10k)
    format (8-column):
    chrom, start, end, ref, depth, ref/alt, base, count
    """
    if isinstance(x, pysam.VariantRecord):
        # ref and alleles
        pos = [x.chrom, x.start, x.stop, x.ref]
        for sample in x.samples:
            xs = x.samples[sample]
            dp = xs.get('DP', None)
            ad = xs.get('AD', None)
            break # only for first sample !!!
        if ad is None:
            msg = ' '.join(list(map(str, pos)))
            print('AD not found in FORMAT: ', msg)
            return None
        for i,j in zip(x.alleles, ad):
            if i in 'ACGT':
                g = 'ref' if i == x.ref else 'alt'
                yield pos + [x.info.get('DP'), g, i, j]


def parse_mpileup(x):
    # check index-ed: .bcf.csi
    if isinstance(x, str):
        if os.path.exists(x) and x.endswith('.bcf'):
            idx = x + '.csi'
            if not os.path.exists(x+'.csi'):
                os.system(f'bcftools index {x}')
            # do
            vcf_in = pysam.VariantFile(x)
            out_txt = os.path.splitext(x)[0]+'.variant.txt'
            with open(out_txt, 'wt') as w:
                for rec in vcf_in.fetch():
                    for s in read_mpileup(rec):
                        w.write('\t'.join(list(map(str, s)))+'\n')
        else:
            print(f'bcf file not exists: {x}')
    else:
        print(f'illegal bcf file, str expect, got {type(x).__name__}')



def main():
    args = vars(get_args().parse_args())
    bam_list = args.pop('bam_list')
    ref = args.pop('ref')
    run_mpileup(bam_list, ref, **args)


if __name__ == '__main__':
    main()


## EOF