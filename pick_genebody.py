#!/usr/bin/env python
#-*- encoding: utf-8 -*-
"""
Deprecated, see `cal_pausing_index.py`

Pick the genebody region, relate to TSS
criteria
- the strongest RNAP2 ChIP-seq signal (TSS region: tss_up + tss_down)
- the distal TSS, if multiple TSShave equal ChIP-seq sigals
- H3K4me3 enrichment (ChIP/input) > 4 (optional) 

input files:
1. TSS records (fixed, strongest TSS), BED6, see pick_tss.py
2. gene records BED6

output files
1. genebody region (BED6)
"""


import os
import sys
import pathlib
import argparse
from xopen import xopen
from hiseq.utils.featurecounts import FeatureCounts
from hiseq.utils.bam import Bam
# from hiseq.utils.utils import Config
# from logging import raiseExceptions


def pick_genebody(gene, tss, **kwargs):
    """
    Parameters:
    -----------
    gene : str
        file, gene records in BED format, original genebody,
        gene-start, gene-end
    tss : str
        file, TSS records in BED format, updated TSS,
        TSS-start, TSS-end
    **kwargs
        required; gb_up, gb_down, tes_extra, out_dir, name, overwrite
    """
    args = {
        'out_dir': None,
        'name': None,
    }
    args.update(kwargs)
    args['out_dir'] = fix_out_dir(args['out_dir'])
    if not isinstance(args['name'], str):
        args['name'] = os.path.basename(os.path.splitext(gene)[0])
    # prepare files
    gene_new = os.path.join(args['out_dir'], args['name']+'.gene.bed')
    # 1. update gene TSS
    update_gene_tss(gene, tss, gene_new)
    # 2. extract genebody file
    gb_file = get_genebody_file(gene_new, fmt='bed', **args)
    # output
    return gb_file


def get_genebody_file(x, fmt='bed', **kwargs):
    """
    x : str
        file, gene record in BED format
    fmt : str
        output format, [bed|gtf], default: [bed]
    **kwargs 
        see arguments in get_tss_region()
        required: out_dir, gb_up, gb_down, tes_extra, name, overwrite
    """
    args = {
        'out_dir': None,
        'gb_up': 250,
        'gb_down': 2250,
        'tes_extra': 0,
        'name': None,
        'overwrite': False,
    }
    args.update(kwargs)
    if not isinstance(args['name'], str):
        args['name'] = os.path.splitext(os.path.basename(x))[0]
    if isinstance(args['gb_down'], int):
        msg2 = '[TSS+{},TSS+{}]'.format(args['gb_up'], args['gb_down'])
        n2 = args['name']+'.gbR.TSS_{}_{}'.format(args['gb_up'], args['gb_down'])
    elif args['gb_down'] == 'TES':
        msg2 = '[TSS+{},TES+{}]'.format(args['gb_up'], args['tes_extra'])
        n2 = args['name']+'.gbR.TSS_{}_TES_{}'.format(args['gb_up'], args['tes_extra'])
    else:
        raiseExceptions('unknown gb_down={}'.format(args['gb_down']))
    # 1. genebody region
    args['out_dir'] = fix_out_dir(args['out_dir'])
    out = os.path.join(args['out_dir'], n2+'.'+fmt)
    print('- get genebody regions: {}'.format(msg2))
    with xopen(out, 'wt') as w:
        for gb_left, gb_right, p in get_genebody_region(x, **args):
            gb = list(map(str, [p[0], gb_left, gb_right]+p[3:6]))
            if fmt == 'gtf':
                gb = bed2gtf(gb, 'gene')
            w.write('\t'.join(gb)+'\n')
    return out


def get_genebody_region(x, gb_up=250, gb_down=2250, tes_extra=0, **kwargs):
    """
    Parameters:
    -----------
    x : str
        file, gene records in BED6 format
    gb_up : int
        to define genebody region, on the right of TSS, default: [250] 
    gb_down : int or str
        as gb_up, on the right of TSS, default: [2250]
        could be 'TES'
    tes_extra : int
        to define the right of genebody region, if gb_down='TES', choose extra 
        region to TES, +/- mean on the right/left of TES, default: [0]

    # depict, how to define the TSS/genebody region
    ## genes on forward strand(+):
                 TSS                                          TES    
    (+)------------|-------------------------------------------|---------->
                   |-----gb_up-----|
                   |---------------------gb_down------|
                                   (.genebody region..)
                   |-----gb_up-----|
                                                  gb_down=TES--|
                                   (..genebody region..........)
                   |-----gb_up-----|
                                                  gb_down=TES--|-tes_extra-|
                                   (..genebody region......................)
    ## genes on reverse strand(-):
                  TES                                        TSS    
    (-)<-----------|------------------------------------------|-----------
                                                |----gb_up----|
                              |----------------------gb_down--|
                              (.genebody region.)
                                                |----gb_up----|
                   |---gb_down=TES
                   (......genebody region.......)
                                                |----gb_up----|
       |-tes_extra-|---gb_down=TES
       (......genebody region...................)
    """
    args = {
        'gb_up': 250,
        'gb_down': 2250,
        'tes_extra': 0,
    }
    args.update(kwargs)
    gb_up = args.get('gb_up')
    gb_down = args.get('gb_down')
    tes_extra = args.get('tes_extra')
    # is_valid_bed(x)
    with xopen(x) as r:
        for l in r:
            p = l.strip().split('\t') # BED6
            s,e = list(map(int, p[1:3])) # start, end
            tss = s + 1 if p[5] == '+' else e
            # 2. genebody region
            gbU = tss + gb_up if p[5] == '+' else tss - gb_up
            if isinstance(gb_down, int):
                gbD = tss + gb_down if p[5] == '+' else tss - gb_down
                gbD = e if gbD > e else s + 1 if gbD < s else gbD # fix border
            elif gb_down == 'TES':
                gbD = e + tes_extra if p[5] == '+' else s - tes_extra + 1
            else:
                raiseExceptions('unknown gb_down={}'.format(gb_down))
            yield (gbU, gbD, p) if p[5] == '+' else (gbD, gbU, p)


def fix_out_dir(x):
    """
    fix out_dir, if not "str", set "cwd()"
    """
    if not isinstance(x, str):
        x = pathlib.Path.cwd()
    x = os.path.abspath(x)
    if not os.path.exists(x):
        os.makedirs(x)
    return x


def bed2gtf(x, feature='gene'):
    """
    Convert BED to GTF format
    x : list
        single bed record, BED6
    """
    if len(x) < 3:
        return None
    n1 = '{}:{}-{}'.format(x[0], x[1], x[2])
    name, strand = [x[3], x[5]] if len(x) > 5 else [n1, '+']
    s, e = x[1:3] # start, end
    s = int(s) + 1 # fix BED-chromstart 0-index
    des = 'gene_id "{}"; gene_name "{}"'.format(name, name)
    gtf = [x[0], 'bed', feature, s, e, '.', strand, '.', des]
    return list(map(str, gtf))


def update_gene_tss(gene, tss, gene_new=None):
    """
    Update the TSS for gene
    """
    # load TSS
    t = {}
    with open(tss) as r:
        for l in r:
            p = l.split('\t')
            t.update({p[3]:p[2]}) # gene_name, TSS
    # for genes
    if not isinstance(gene_new, str):
        gene_new = os.path.splitext(gene)[0]+'.TSS_fixed.bed'
    with open(gene) as r, open(gene_new, 'wt') as w:
        for l in r:
            p = l.strip().split('\t')
            b = p.copy()
            s0 = int(p[1])+1 if p[5] == '+' else int(p[2]) # current
            s1 = t.get(p[3], s0) # 1-indexed
            s1 = int(s1)
            b[1], b[2] = (s1-1, p[2]) if p[5] == '+' else (p[1], s1)
            b = list(map(str, b))
            w.write('\t'.join(b)+'\n')
    return gene_new


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tss-bed', dest='tss_bed', required=True, 
        help='The TSS record in BED6 format, see "format_tss.py" for help')
    parser.add_argument('-g', '--gene-bed', dest='gene_bed', required=True,
        help='The genebody in BED format,')
    parser.add_argument('-o', dest='out_dir', required=False, 
        help='The out_dir')
    parser.add_argument('-U', '--gb-up', dest='gb_up', type=int,
        default=250, help='for genebody region, downstream of TSS, default: [250]')
    parser.add_argument('-D', '--gb-down', dest='gb_down', type=int,
        default=2250, help='for genebody region, downstream of TSS, default: [2250]')
    parser.add_argument('-X', '--tes-extra', dest='tes_extra', type=int,
        default=0, help='for genebody region, downstream of TES, default: [0]')
    parser.add_argument('-n', '--name', dest='name', type=str,
        default=None, help='name of the output files, default: [auto]')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite exists file')
    return parser


def main():
    args = vars(get_args().parse_args())
    gene = args.pop('gene_bed')
    tss = args.pop('tss_bed')
    pick_genebody(gene, tss, **args)


if __name__ == '__main__':
    main()

# EOF