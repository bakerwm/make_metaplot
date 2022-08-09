#!/usr/bin/env python
#-*- encoding: utf-8 -*-
"""
Pick TSS and genebody regions

option-1: choose gene-start as "TSS"
arguments: gene_bed, 

option-2: pick the "TSS" with strongest Pol II ChIP signal
arguments: gene_bed, tss_bed, bam, ...

common arguments:
out_dir, tss_up, tss_down, gb_up, gb_down, tes_extra, name, overwrite


## bugs:
to-do: 
1. remove duplicate genes in GTF annotation (Ensembl, release-102)
"""


import os
import sys
import re
import pathlib
import argparse
from xopen import xopen
from hiseq.utils.featurecounts import FeatureCounts
from hiseq.utils.bam import Bam
from hiseq.utils.file import symlink_file
from hiseq.utils.utils import Config
from logging import raiseExceptions


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
    args_fc['out_dir'] = os.path.join(args_fc['out_dir'], 'count_regsions')
    # print('!A-2', tss_r, bam, args_fc)
    tss_fc = count_region(tss_r, bam, **args_fc)
    gb_fc = count_region(gb_r, bam, **args_fc)
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


################################################################################
# option-1: pick TSSs with the strongest Pol II ChIP sigal
# required: gene_bed, tss_bed, tss_bam
def pick_tss_genebody(gene_bed, tss_bed=None, tss_bam=None, **kwargs):
    args = {
        'out_dir': None
    }
    args.update(kwargs)
    overwrite = args.get('overwrite', False)
    args['out_dir'] = fix_out_dir(args['out_dir'])
    if tss_bed is not None and tss_bam is not None:
        # 1. pick the strongest TSS
        args1 = args.copy()
        args1['out_dir'] = os.path.join(args1['out_dir'], 'top_TSS') # subdir
        fix_out_dir(args1['out_dir'])
        tss_top = pick_top_tss(tss_bed, tss_bam, **args1)
        # 2. update TSS for gene_bed
        args2 = args.copy()
        args2['out_dir'] = os.path.join(args2['out_dir'], 'updated_gene')
        fix_out_dir(args2['out_dir'])
        gene_new = os.path.join(args2['out_dir'], os.path.basename(gene_bed))
        gene_new = update_tss(gene_bed, tss_top, gene_new, overwrite)
    else:
        gene_new = gene_bed
    # 3. pick tss
    tss_f = get_tss_file(gene_new, **kwargs)
    # 4. pick genebody
    gb_f = get_genebody_file(gene_new, **kwargs)
    return [tss_f, gb_f]


################################################################################
def get_tss_file(x, fmt='bed', **kwargs):
    """
    x : str
        gene record in BED format, file
    fmt : str
        output format, [bed|gtf], default: [bed]
    **kwargs 
        see arguments in get_tss_region()
        required: out_dir, tss_up, tss_down, name, overwrite
    """
    args = {
        'out_dir': None,
        'tss_up': 150,
        'tss_down': 150,
        'name': None,
        'overwrite': False,
    }
    args.update(kwargs)
    if not isinstance(args['name'], str):
        args['name'] = os.path.basename(os.path.splitext(x)[0])
    msg1 = '[TSS-{},TSS+{}]'.format(args['tss_up'], args['tss_down'])
    n1 = args['name']+'.TSSR.{}_TSS_{}'.format(args['tss_up'], args['tss_down'])
    # 1. TSS region
    args['out_dir'] = fix_out_dir(args['out_dir'])
    out = os.path.join(args['out_dir'], n1+'.'+fmt)
    print('- get TSS regions: {}'.format(msg1))
    with xopen(out, 'wt') as w:
        # for tss_left, tss_right, p in get_tss_region(x, **args):
        #     tss = list(map(str, [p[0], tss_left, tss_right]+p[3:6]))
        for b in get_tss_region(x, **args):
            if fmt == 'gtf':
                b = bed2gtf(b, 'gene')
            w.write('\t'.join(b)+'\n')
    return out


def get_tss_region(x, tss_up=150, tss_down=150, **kwargs):
    """
    Parameters:
    -----------
    x : str
        gene records in BED6 format; gene record
    tss_up : int
        upstream border of TSS region, distance to TSS site, default: [150]
    tss_down : int
        downstream border of TSS region, distance to TSS site, default: [150]
    # depict, how to define the TSS region
    # genes on forward strand(+):
                 TSS                                          TES    
    (+)------------|----------------[gene]--------------------|--------->
        |--tss_up--|--tss_down--|
        (.......TSS rgion.......)

    ## genes on reverse strand(-):
                  TES                                        TSS    
    (-)<-----------|----------------[gene]--------------------|-----------
                                                 |--tss_down--|--tss_up--|
                                                 (.......TSS rgion.......)
    """
    print('!B-2', kwargs, type(tss_up), type(tss_down))
    # is_valid_bed(x)
    with xopen(x) as r:
        for l in r:
            p = l.strip().split('\t') # BED6
            s,e = list(map(int, p[1:3])) # start, end
            s,e = (s-tss_up, s+tss_down+1) if p[5] == '+' else (e-tss_down-1, e+tss_up)
            if s < 0:
                s = 0
            b = [p[0], s, e] + p[3:6]
            yield list(map(str, b)) # return BED record


################################################################################
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
    # print('!A-3', args['gb_down'], type(args['gb_down']))
    # sys.exit(1)
    if not isinstance(args['name'], str):
        args['name'] = os.path.splitext(os.path.basename(x))[0]
    if isinstance(args['gb_down'], int):
        msg2 = '[TSS+{},TSS+{}]'.format(args['gb_up'], args['gb_down'])
        n2 = args['name']+'.gbR.TSS_{}_{}'.format(args['gb_up'], args['gb_down'])
        # print('!A-4', n2, msg2)
    elif args['gb_down'] == 'TES':
        msg2 = '[TSS+{},TES+{}]'.format(args['gb_up'], args['tes_extra'])
        n2 = args['name']+'.gbR.TSS_{}_TES_{}'.format(args['gb_up'], args['tes_extra'])
        # print('!A-5', n2, msg2)
    else:
        raiseExceptions('unknown gb_down={}'.format(args['gb_down']))
        # print('!A-6', args['gb_down'])
    # 1. genebody region
    args['out_dir'] = fix_out_dir(args['out_dir'])
    out = os.path.join(args['out_dir'], n2+'.'+fmt)
    print('- get genebody regions: {}'.format(msg2))
    with xopen(out, 'wt') as w:
        # for gb_left, gb_right, p in get_genebody_region(x, **args):
        #     gb = list(map(str, [p[0], gb_left, gb_right]+p[3:6]))
        # print('!B1, {}'.format(args.get('gb_down', None)))
        for b in get_genebody_region(x, **args):
            if fmt == 'gtf':
                b = bed2gtf(b, 'gene')
            w.write('\t'.join(b)+'\n')
    return out


def get_genebody_region(x, **kwargs):
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
            else:
                raiseExceptions('unknown gb_down={}'.format(gb_down))
            if s > e:
                print('line-{}, illegal BED: {}'.format(i, l.strip()))
            b = p.copy()
            b[1], b[2] = list(map(str, [s, e])) #update
            yield b


################################################################################
# choose the strongest TSS
def pick_top_tss(x, bam, **kwargs):
    """
    x : str
        TSS records in BED6 format, see script "format_tss.py" for help
    bam : str
        The RNAP2 ChIP alignment file, to determine TSS site
    out_dir : str
        Directory to save final records
    # y : str
    #     The output TSS, only one for each gene
    tss_up : int
        TSS region, upstream of TSS, default: [150]
    tss_down : int
        TSS region, downstream of TSS, default: [150]
    ...
    """
    print('- Pick the strongest TSS')
    args = {
        'name': None,
    }
    args.update(kwargs)
    overwrite = args.get('overwrite', False)
    # 1. extract TSS region
    tss_f = get_tss_file(x, fmt='bed', **args)
    # 2. count reads on TSS regions
    args['name'] = os.path.basename(os.path.splitext(tss_f)[0]) # update name
    tss_fc = count_region(x=tss_f, bam=bam, **args)
    # 3. pick TSS sites
    tss_top = pick_distal_tss(tss_fc, overwrite)
    # 4. convert TSS to single-base BED format
    tss_top_fixed = recover_tss(tss_top) # retrieve from total TSS, by gene_name
    # 5. save top_tss to main directory
    tss = os.path.join(args['out_dir'], os.path.basename(tss_top_fixed))
    symlink_file(tss_top_fixed, tss)
    return tss #_top_fixed


def pick_distal_tss(x, overwrite=False):
    """
    Pick the distal TSS for genes with multiple TSSs that equal in RNAP2 signal
    x : str
        otuput of featureCounts, read counts in column-7
    """
    x_new = os.path.splitext(x)[0]+'.top_TSS.bed'
    if os.path.exists(x_new) and overwrite is False:
        print('pick_distal_tss() skipped, file exists: {}'.format(x_new))
        return x_new
    d = {}
    for p in load_fc(x):
        # dup gene names in multiple chromosome
        g,c,_ = p[0].split(':', 2) # fix gene_name:tss
        k = '{}:{}'.format(g,c) # gene+chr
        dg = d.get(k, [])
        dg.append(p)
        d.update({k:dg})
    # pick strongest tss
    print('Number of TSSs for gene: {}'.format(len(d)))
    # choose distal tss, for multiple TSSs with equal counts
    with open(x_new, 'wt') as w:
        for k,fc in d.items(): # k=gene+chr
            g,c = k.split(':', 1)
            n = [float(i[6]) for i in fc] # count in column-7
            max_idx = [i for i,j in enumerate(n) if j == max(n)] # max count
            if len(max_idx) == 1:
                max_fc = fc[max_idx[0]]
            else:
                fc2 = [j for i,j in enumerate(fc) if i in max_idx]
                max_fc = pick_distal_tss2(fc2)
            # fix BED start, 0-index
            max_fc[2] = str(int(max_fc[2])-1) 
            t = max_fc[1:4]+[g, 254, max_fc[4]] # BED format
            t = list(map(str, t))
            w.write('\t'.join(t)+'\n')
    return x_new


def pick_distal_tss2(x):
    """
    x : list
        list of featureCounts records, list of list
    determine the distal TSS
    """
    try:
        s = [int(i[2]) for i in x] # gene_id, chr, start, end, strand, length, count
    except:
        print('!A-3', x)
        sys.exit(1)
    distal_s = min(s) if x[0][4] == '+' else max(s)
    return x[s.index(distal_s)]


def recover_tss(x, overwrite=False):
    """
    Recover the TSS region to TSS site
    parse tss_up and tss_down from filename: ...TSSR.150_TSS_150.top_TSS
    Parameters:
    -----------
    x : str
        file, fixed TSS in BED format, see pick_distal_tss()
        suffix: ".TSSR.150_TSS_150.top_TSS"

    # depict, how to define the TSS region
    # genes on forward strand(+):
                 TSS                                          TES    
    (+)------------|----------------[gene]--------------------|--------->
        |--tss_up--|--tss_down--|
        (.......TSS rgion.......)

    ## genes on reverse strand(-):
                  TES                                        TSS    
    (-)<-----------|----------------[gene]--------------------|-----------
                                                 |--tss_down--|--tss_up--|
                                                 (.......TSS rgion.......)
    write the final results to file: ".TSS.bed"
    """
    x_new = re.sub('.TSSR.\d+_TSS_\d+\.top_TSS', '.TSS', x, flags=re.IGNORECASE)
    if os.path.exists(x_new) and overwrite is False:
        print('recover_tss() skipped, file exists: {}'.format(x_new))
        return x_new
    # if not os.path.exists(x):
    #     print('file not exists: {}'.format(x))
    #     sys.exit(1)
    # fix filename, remove suffix '.TSSR.150_TSS_150'
    g = re.search('\.TSSR\.(\d+)_TSS_(\d+)', x)
    if g is None:
        print('x, expect [.TSSR.150_TSS_150], got {}'.format(x))
        sys.exit(1)
    tss_up, tss_down = (int(g.group(1)), int(g.group(2)))
    # new name
    i = 0
    with open(x_new, 'wt') as w, open(x) as r:
        for l in r:
            i += 1
            p = l.strip().split('\t')
            s,e = list(map(int, p[1:3])) # start, end
            s,e = (s+tss_up, e-tss_down) if p[5] == '+' else (s+tss_down, e-tss_up)
            if s >= e:
                s = e - 1 # one the edges
                print('line-{}, illegal BED, {}, {}'.format(i, s, e))
            p[1], p[2] = list(map(str, [s, e]))
            w.write('\t'.join(p)+'\n')
    return x_new


# count reads on TSS region
def count_region(x, bam, **kwargs):
    """
    Count reads on genebody region using featureCounts
    x : str
        gene record in BED/GTF format, file
    kwargs : optional arguments
        required, out_dir, name, overwrite
    """
    args = {
        'out_dir': None,
        'name': None,
        'overwrite': False,
    }
    args.update(kwargs)
    args['out_dir'] = fix_out_dir(args['out_dir'])
    bname = os.path.basename(os.path.splitext(bam)[0]) # bam name
    args['out_dir'] = os.path.join(args['out_dir'], bname)
    if not isinstance(args['name'], str):
        args['name'] = os.path.splitext(os.path.basename(x))[0]
    # 1. get gtf file
    if x.endswith('.gtf'):
        gtf = x
    elif x.endswith('.bed'):
        gtf = x.replace('.bed', '.gtf')
        with open(x) as r, open(gtf, 'wt') as w:
            for l in r:
                p = l.strip().split('\t')
                g = bed2gtf(p, feature='gene')
                w.write('\t'.join(g)+'\n')
    else:
        sys.exit('unknown x={}'.format(x))
    # 2. count genebody region
    args_fc = {
        'gtf': gtf,
        'bam_list': bam,
        'outdir': args['out_dir'],
        'prefix': args['name']+'.txt',
        'feature_type': 'gene',
    }
    # index bam file
    Bam(bam).index()
    fc = FeatureCounts(**args_fc)
    fc.run()
    return fc.count_txt


def load_fc(x):
    with xopen(x) as r:
        for l in r:
            if l.startswith('#') or l.startswith('Geneid'):
                continue
            p = l.strip().split('\t')
            yield p


################################################################################
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


def update_tss(in_bed, tss_bed, gene_new=None, overwrite=False):
    """
    Update the TSS for BED
    Parameters:
    -----------
    in_bed : str
        file, gene records in BED6 format
    tss_bed : str
        file, TSS records in BED6 format, single-base
    gene_new : str, None
        str or None, save the updated gene-records
    """
    print('- Update TSS for genes')
    # for genes
    if not isinstance(gene_new, str):
        gene_new = os.path.splitext(in_bed)[0]+'.TSS_updated.bed'
    if os.path.exists(gene_new) and overwrite is False:
        print('update_tss() skipped, file exists: {}'.format(gene_new))
        return gene_new
    # load TSS
    t = {}
    i = 0
    with open(tss_bed) as r:
        for l in r:
            i += 1
            p = l.split('\t')
            s,e = list(map(int, p[1:3])) # start, end
            if e - s > 1:
                raiseExceptions('line-{}, not a single-base, {}'.format(i, l))
            t.update({p[3]:p[2]}) # gene_name, TSS
    # iterate genes
    with open(in_bed) as r, open(gene_new, 'wt') as w:
        for l in r:
            p = l.strip().split('\t')
            b = p.copy()
            s, e = list(map(int, b[1:3])) # start, end
            if s >= e: # !!! illegal genes, with same id !!! 
                print('illegal gene coordinates, {}'.format(p))
                continue
            tss = s + 1 if b[5] == '+' else e # current TSS
            tss_new = t.get(p[3], tss) # updated TSS
            tss_new = int(tss_new)
            b[1], b[2] = [tss_new - 1, e] if p[5] == '+' else [s, tss_new]
            if b[1] >= b[2]: # !!! illegal genes, with same id
                print('illegal gene coordinates, {}'.format(p))
                continue
            b = list(map(str, b))
            w.write('\t'.join(b)+'\n')
    return gene_new


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gene-bed', dest='gene_bed', required=True,
        help='The genebody in BED format,')
    parser.add_argument('-t', '--tss-bed', dest='tss_bed', required=False, 
        help='The TSS record in BED6 format, see "format_tss.py" for help')
    parser.add_argument('-b', '--tss-bam', dest='tss_bam', required=False,
        help='The alignemnt file of Pol2 ChIP, in BAM format')
    parser.add_argument('-o', dest='out_dir', required=False, 
        help='The out_dir')
    parser.add_argument('-u', '--tss-up', dest='tss_up', type=int,
        default=150, help='for TSS region, upstream of TSS, default: [150]')
    parser.add_argument('-d', '--tss-down', dest='tss_down', type=int,
        default=150, help='for TSS region, downstream of TSS, default: [150]')
    parser.add_argument('-U', '--gb-up', dest='gb_up', type=int,
        default=250, help='for genebody region, downstream of TSS, default: [250]')
    parser.add_argument('-D', '--gb-down', dest='gb_down', 
        default=2250, help='for genebody region, downstream of TSS, default: [2250]')
    parser.add_argument('-X', '--tes-extra', dest='tes_extra', type=int,
        default=0, help='for genebody region, downstream of TES, default: [0]')
    parser.add_argument('-n', '--name', dest='name', type=str,
        default=None, help='name of the output files, default: [auto]')
    parser.add_argument('-x', '--pi-ver', dest='pi_ver', type=int, default=1,
        help='version of pausing index, 1=rpk, 2=rpbm, 3=rpm, default: [1]')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite exists file')
    return parser


def main():
    args = vars(get_args().parse_args())
    # fix gb_down
    if re.match('^[0-9]+$', args['gb_down']):
        args['gb_down'] = int(args['gb_down']) # convert to int
    gene = args.pop('gene_bed')
    tss_bed = args.pop('tss_bed', None)
    tss_bam = args.pop('tss_bam', None)
    tss_r, gb_r = pick_tss_genebody(gene, tss_bed, tss_bam, **args)
    cal_pausing_index(tss_r, gb_r, tss_bam, **args)


if __name__ == '__main__':
    main()

# EOF