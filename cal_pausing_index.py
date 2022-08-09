#!/usr/bin/env python
#-*- encoding: utf-8 -*-
"""
Calculate Pausing index:

TSS region: -150bp to +150bp (arguments)
genebody region: +150bp to 2250bp (arguments)

# remove short genes

common arguments:
out_dir, tss_up, tss_down, gb_up, gb_down, tes_extra, name, overwrite
"""


import os
import sys
import re
import pathlib
import argparse
from xopen import xopen
from hiseq.utils.featurecounts import FeatureCounts
from hiseq.utils.bam import Bam
from hiseq.utils.utils import Config
from logging import raiseExceptions
from utils import fix_out_dir, log, symlink_file
from pick_tss import get_tss_file, count_region, load_fc


def cal_pausing_index(tss_r, gb_r, bam, **kwargs):
    """
    Calculate the Pausing index (PI):
    p_idx = tss_density / genebody_density

    version-1: 2015 mol cell, Karen Adelman lab
    score = rpk (reads/kb or read_pairs/kb)
    version-2: 2016 genome biol, William lab 
    score = rpbm (rpm / bp)
    version-3: 2018 cell, Ali Shilatifard lab
    score = rpm
    
    required parameters
    see get_args()
    bed, bam, out_dir, tss_up, tss_down, gb_up, gb_down, tes_extra, name
    """
    args = {
        'scale': 1,
        'pi_ver': 1,
        'overwrite': False
    }
    args.update(kwargs)
    # 1. count reads on TSS/genebody
    args_fc = args.copy()
    # print('!A-2', tss_r, bam, args_fc)
    args_fc['prefix'] = None
    tss_fc = count_region(tss_r, bam, **args_fc)
    try:
        gb_fc = count_region(gb_r, bam, **args_fc)
    except:
        print('!B-1', args_fc)
        sys.exit(1)
    # 2. calculate pausing index
    args['scale'] = cal_norm_scale(bam, os.path.dirname(tss_fc), force=False)
    pidx_f = cal_pidx(tss_fc, gb_fc, **args)
    print('save pausing index: {}'.format(pidx_f))
    return pidx_f


def cal_pidx(tss_count, gb_count, scale=1, pi_ver=1, overwrite=False, **kwargs):
    """
    Calculate Pausing index
    pidx = tss_density / genebody_density

    Parameters:
    -----------
    tss_count : str
        file, featureCounts output of TSS regions
    gb_count : str
        file, featureCounts output of genebody regions
    scale : float
        the scale to normalize the density, default: [1]
    pi_ver : int
        how to calculate the p_idx, [1,2,3], default: [1]
        1="rpk, reads/kilo base", 2="rpbm, rpm/base", 
        3="rpm, reads per million"
    overwrite : bool
        overwrite the exists files
    """
    # pidx version
    s = 'rpk' if pi_ver==1 else 'rpbm' if pi_ver==2 else 'rpm'
    # output file
    name = 'pausing_index.{}.txt'.format(s)
    out_f = os.path.join(os.path.dirname(tss_count), name)
    # save to file
    if False: # os.path.exists(out_f) and overwrite is False:
        pass
    else:
        # load densities
        tss_density = load_density(tss_count, scale, pi_ver)
        gb_density = load_density(gb_count, scale, pi_ver)
        with open(out_f, 'wt') as w:
            for k,v_tss in tss_density.items():
                v_gb = gb_density.get(k, [1])
                try:
                    p = '{:.4f}'.format(float(v_tss[-1]) / float(v_gb[-1]))
                except ZeroDivisionError as e:
                    print('gene:{}, divided by zero'.format(k))
                    p = 0
                # print('!A-1', p, v_tss[-1], v_gb[-1])
                # output
                f = v_tss[:8] + v_gb[1:8] + [p]
                f = list(map(str, f))
                w.write('\t'.join(f)+'\n')
    return out_f


def load_density(x, scale=1, pi_ver=1):
    """
    load density from featureCounts output
    v1: reads/kb
    v2: rpm/bp
    v3: rpm
    input: output of featureCounts.txt
    Geneid Chr Start End Strand Length bam [bam ...]
    column-1: Geneid
    column-2: Chr
    column-3: Start
    column-4: End
    column-5: Strand
    column-6: Length
    column-7: bam    
    output: density
    """
    cal_density = cal_density_v1 if pi_ver == 1 else cal_density_v2 if pi_ver == 2 else cal_density_v3
    d = {}
    with xopen(x) as r:
        for l in r:
            if l.startswith('#') or l.startswith('Geneid'):
                continue
            p = l.strip().split('\t')
            n = cal_density(p[6], p[5], scale)
            d.update({
                p[0]: p[:7] + [n],
            })
    return d


def cal_density_v1(n_reads, gene_size, scale=1):
    # version-1: reads/kb
    return '{:.6f}'.format(float(n_reads) / (int(gene_size) / 1000)) 


def cal_density_v2(n_reads, gene_size, scale=1):
    # version-2: rpm/bp 
    return '{:.6f}'.format((float(n_reads) * scale) / int(gene_size)) 


def cal_density_v3(n_reads, gene_size, scale=1):
    # version-3: rpm
    try:
        d = '{:.6f}'.format((float(n_reads) * scale))
    except:
        print('!A-3', n_reads, gene_size, scale)
        d = 0
    return d


def cal_norm_scale(bam, out_dir=None, norm=1000000, force=False):
    """
    Calculate the norm scale for bam file
    output: out_dir/scale.json
    format:
    {
        name: bam_name,
        total: 8000000,
        norm: 1000000,
        scale: 0.1250,
    }
    """
    print('- cal bam norm scale: {}'.format(os.path.basename(bam)))
    # prepare output file
    # if isinstance(out_dir, str):
    #     out_dir = fix_out_dir(out_dir)
    #     out_json = os.path.join(out_dir, 'scale.json')
    # else:
    #     out_json = None        
    out_dir = fix_out_dir(out_dir)
    out_json = os.path.join(out_dir, 'scale.json')
    # load norm scale
    if os.path.exists(out_json) and force is False:
        d = Config().load(out_json)
        s = d.get('scale', 1) # default
    else:
        m = Bam(bam).count(reads=False) # read_pairs, or se
        s = '{:.4f}'.format(norm / m)
        d = {
            'name': os.path.basename(os.path.splitext(bam)[0]),
            'total': m,
            'norm': norm,
            'scale': float(s),
        }
        if out_json is not None:
            Config().dump(d, out_json)
    return s


def get_genebody_file(x, fmt='bed', **kwargs):
    """
    Parameters:
    -----------
    x : str
        file, gene record in BED format
    fmt : str
        output format, [bed|gtf], default: [bed]
    **kwargs 
        see arguments in get_tss_region()
        required: out_dir, gb_up, gb_down, tes_extra, prefix, overwrite
    """
    args = {
        'out_dir': None,
        'gb_up': 250,
        'gb_down': 2250,
        'tes_extra': 0,
        'prefix': None,
        'overwrite': False,
    }
    args.update(kwargs)
    # print('!A-3', args['gb_down'], type(args['gb_down']))
    # sys.exit(1)
    # update files
    args['out_dir'] = fix_out_dir(args['out_dir'])
    if not isinstance(args['prefix'], str):
        args['prefix'] = os.path.splitext(os.path.basename(x))[0]
    if isinstance(args['gb_down'], int):
        msg = 'genebody region: --TSS--[({}bp)--({}bp)]--'.format(args['gb_up'], args['gb_down'])
        n2 = args['prefix']+'.gbR.TSS_{}_{}'.format(args['gb_up'], args['gb_down'])
        # print('!A-4', n2, msg)
    elif isinstance(args['gb_down'], float):
        msg = 'genebody region: --TSS--[({}%)--({}%)]--'.format(args['gb_up']*100, args['gb_down']*100)
        n2 = args['prefix']+'.gbR.TSS_{}_{}'.format(args['gb_up'], args['gb_down'])
        # print('!A-4', n2, msg)
    elif args['gb_down'] == 'TES':
        msg = 'genebody region: --TSS--[({}bp)--TES--({}bp)]--'.format(args['gb_up'], args['tes_extra'])
        n2 = args['prefix']+'.gbR.TSS_{}_TES_{}'.format(args['gb_up'], args['tes_extra'])
        # print('!A-5', n2, msg)
    else:
        raiseExceptions('unknown gb_down={}'.format(args['gb_down']))
        # print('!A-6', args['gb_down'])
    # 1. genebody region
    out = os.path.join(args['out_dir'], n2+'.'+fmt)
    print(msg)
    if os.path.exists(out) and args.get('overwrite', False) is False:
        print('file exists: {}'.format(out))
    with xopen(out, 'wt') as w:
        for b in get_genebody_region(x, **args):
            if fmt == 'gtf':
                b = bed2gtf(b, 'gene')
            w.write('\t'.join(b[:6])+'\n')
    return out

def get_genebody_region(x, **kwargs):
    """
    Extract genebody region by ratio/bins
    gb_up/gb_down is float [0-1]

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
    i = 0
    with xopen(x) as r:
        for l in r:
            i += 0
            p = l.strip().split('\t') # BED6
            s,e = list(map(int, p[1:3])) # start, end
            if isinstance(gb_down, int):
                s,e = (s+gb_up, s+gb_down) if p[5] == '+' else (e-gb_down-1, e-gb_up)
            elif gb_down == 'TES':
                s,e = (s+gb_up, e+tes_extra) if p[5] == '+' else (s-tes_extra, e-gb_up)
            elif isinstance(gb_down, float): # [0-1]
                if p[5] == '+':
                    e = int(s + (e - s ) * gb_down)
                else:
                    s = int(e - (e - s) * gb_down)
            elif isinstance(gb_up, float):
                if p[5] == '+':
                    s = int(s + (e - s) * gb_up)
                else:
                    e = int(e - (e - s) * gb_dup)
            else:
                raiseExceptions('unknown gb_down={}'.format(gb_down))
            if s > e:
                print('line-{}, illegal BED: {}'.format(i, l.strip()))
            b = p.copy()
            # skip short genes
            if s >= e:
                continue # !!!
            b[1], b[2] = list(map(str, [s, e])) #update

            yield b


# def get_genebody_region(x, **kwargs):
#     """
#     Parameters:
#     -----------
#     x : str
#         file, gene records in BED6 format
#     gb_up : int
#         to define genebody region, on the right of TSS, default: [250] 
#     gb_down : int or str
#         as gb_up, on the right of TSS, default: [2250]
#         could be 'TES'
#     tes_extra : int
#         to define the right of genebody region, if gb_down='TES', choose extra 
#         region to TES, +/- mean on the right/left of TES, default: [0]

#     # depict, how to define the TSS/genebody region
#     ## genes on forward strand(+):
#                  TSS                                          TES    
#     (+)------------|-------------------------------------------|---------->
#                    |-----gb_up-----|
#                    |---------------------gb_down------|
#                                    (.genebody region..)
#                    |-----gb_up-----|
#                                                   gb_down=TES--|
#                                    (..genebody region..........)
#                    |-----gb_up-----|
#                                                   gb_down=TES--|-tes_extra-|
#                                    (..genebody region......................)
#     ## genes on reverse strand(-):
#                   TES                                        TSS    
#     (-)<-----------|------------------------------------------|-----------
#                                                 |----gb_up----|
#                               |----------------------gb_down--|
#                               (.genebody region.)
#                                                 |----gb_up----|
#                    |---gb_down=TES
#                    (......genebody region.......)
#                                                 |----gb_up----|
#        |-tes_extra-|---gb_down=TES
#        (......genebody region...................)
#     """
#     args = {
#         'gb_up': 250,
#         'gb_down': 2250,
#         'tes_extra': 0,
#     }
#     args.update(kwargs)
#     gb_up = args.get('gb_up')
#     gb_down = args.get('gb_down')
#     tes_extra = args.get('tes_extra')
#     # is_valid_bed(x)
#     i = 0
#     with xopen(x) as r:
#         for l in r:
#             i += 0
#             p = l.strip().split('\t') # BED6
#             s,e = list(map(int, p[1:3])) # start, end
#             if isinstance(gb_down, int):
#                 s,e = (s+gb_up, s+gb_down) if p[5] == '+' else (e-gb_down-1, e-gb_up)
#             elif gb_down == 'TES':
#                 s,e = (s+gb_up, e+tes_extra) if p[5] == '+' else (s-tes_extra, e-gb_up)
#             else:
#                 raiseExceptions('unknown gb_down={}'.format(gb_down))
#             if s > e:
#                 print('line-{}, illegal BED: {}'.format(i, l.strip()))
#             b = p.copy()
#             # skip short genes
#             if s >= e:
#                 continue # !!!
#             b[1], b[2] = list(map(str, [s, e])) #update

#             yield b


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


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gene-bed', dest='gene_bed', required=True,
        help='The genebody in BED format,')
    parser.add_argument('-b', '--bam', dest='bam', required=False,
        help='The alignemnt file of Pol2 ChIP, in BAM format')
    parser.add_argument('-o', dest='out_dir', required=False, 
        help='The out_dir')
    parser.add_argument('-u', '--tss-up', dest='tss_up', type=int,
        default=150, help='for TSS region, upstream of TSS, default: [150]')
    parser.add_argument('-d', '--tss-down', dest='tss_down', type=int,
        default=150, help='for TSS region, downstream of TSS, default: [150]')
    parser.add_argument('-U', '--gb-up', dest='gb_up', type=str,
        default=250, help='for genebody region, downstream of TSS, default: [250]')
    parser.add_argument('-D', '--gb-down', dest='gb_down', type=str,
        default='2250', help='for genebody region, downstream of TSS, default: [2250]')
    parser.add_argument('-X', '--tes-extra', dest='tes_extra', type=int,
        default=0, help='for genebody region, downstream of TES, default: [0]')
    parser.add_argument('-n', '--prefix', dest='prefix', type=str,
        default=None, help='name of the output files, default: [auto]')
    parser.add_argument('-x', '--pi-ver', dest='pi_ver', type=int, default=1,
        help='version of pausing index, 1=rpk, 2=rpbm, 3=rpm, default: [1]')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite exists file')
    return parser


def main():
    args = vars(get_args().parse_args())
    if re.match('^[0-9.]+$', args['gb_down']):
        args['gb_down'] = eval(args['gb_down']) # to int or float
    if re.match('^[0-9.]+$', args['gb_up']):
        args['gb_up'] = eval(args['gb_up']) # to int or float
        # args['gb_down'] = int(args['gb_down']) # convert to int
    gene = args.pop('gene_bed')
    bam = args.pop('bam', None)
    args['out_dir'] = fix_out_dir(args['out_dir'])
    # print('!A-3', args)
    # 1. region dir
    args1 = args.copy()
    args1.update({'out_dir': os.path.join(args['out_dir'], 'region_files')})
    tss_f = get_tss_file(gene, **args1)
    gb_f = get_genebody_file(gene, **args1)
    # 2. count dir
    args2 = args.copy()
    args2.update({'out_dir': os.path.join(args['out_dir'], 'count_files')})
    cal_pausing_index(tss_f, gb_f, bam, **args2)


if __name__ == '__main__':
    main()

# EOF