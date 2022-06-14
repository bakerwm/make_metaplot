#!/usr/bin/env python
"""
extract the TSS regions and genebody regions
Parameters:
TSS: tss_up=150, tss_down=150
genebody: gb_up = 250 (relate to TSS), gb_down = 2250/TES, tes_extra=0 (eg:TES+3k)

To-Do
1. choose the correct TSS site, for genes with multiple TSSs.
"""

from logging import raiseExceptions
import os
import sys
import pathlib
import argparse
from xopen import xopen
from hiseq.utils.featurecounts import FeatureCounts
from hiseq.utils.bam import Bam
from hiseq.utils.utils import Config


def cal_pausing_index(**kwargs):
    """
    Calculate the Pausing index (PI):
    PI = (TSS region) / (genebody region)
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
    args = kwargs.copy()
    # 1. count reads on TSS/genebody
    fc_tss = count_region(x=args['bed'], region='TSS', **args)
    fc_gb = count_region(x=args['bed'], region='genebody', **args)
    # 2. calculate TSS / genebody
    overwrite = args.get('overwrite', False)
    pi_ver = args.get('pi_ver', 1) # norm version [1-3]
    scale = cal_norm_scale(args['bam'], os.path.dirname(fc_tss), force=False)
    p_idx_f = cal_pidx(fc_tss, fc_gb, scale, pi_ver, overwrite)
    print('pausing index: {}'.format(p_idx_f))


def cal_pidx(fc_tss, fc_gb, scale, pi_ver, overwrite):
    """
    Calculate Pausing index
    """
    # output suffix
    s = 'rpk' if pi_ver==1 else 'rpbm' if pi_ver==2 else 'rpm'
    p_idx = 'pausing_index.{}.txt'.format(s)
    p_idx_f = os.path.join(os.path.dirname(fc_tss), p_idx)
    # save to file
    if os.path.exists(p_idx_f) and overwrite is False:
        pass
        # print('file exists: {}'.format(p_idx_f))
    else:
        tss_density = load_density(fc_tss, scale, pi_ver)
        gb_density = load_density(fc_gb, scale, pi_ver)
        with open(p_idx_f, 'wt') as w:
            for k,v_tss in tss_density.items():
                v_gb = gb_density.get(k, 1)
                try:
                    p = '{:.4f}'.format(float(v_tss[-1]) / float(v_gb[-1]))
                except ZeroDivisionError as e:
                    print('gene:{}, divided by zero'.format(k))
                    p = 0
                # output
                f = v_tss[:8] + v_gb[1:8] + [p]
                f = list(map(str, f))
                w.write('\t'.join(f)+'\n')
    return p_idx_f


def load_density(x, scale=1, pi_ver=1):
    """
    Calculate the region density for TSS, genebody
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
    return '{:.6f}'.format((float(n_reads) * scale)) 


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
    out_dir = fix_out_dir(out_dir)
    out_json = os.path.join(out_dir, 'scale.json')
    if os.path.exists(out_json) and force is False:
        d = Config().load(out_json)
        s = d.get('scale', 1) # default
    else:
        m = Bam(bam).count(reads=False) # read_pairs, or se
        s = '{:.4f}'.format(norm / m)
        s = float(s)
        d = {
            'name': os.path.basename(os.path.splitext(bam)[0]),
            'total': m,
            'norm': norm,
            'scale': s
        }
        Config().dump(d, out_json)
    return s


def get_tss_region(x, tss_up=150, tss_down=150, **kwargs):
    """
    Parameters:
    -----------
    x : str
        gene records in BED6 format
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
    # is_valid_bed(x)
    with xopen(x) as r:
        for l in r:
            p = l.strip().split('\t') # BED6
            s, e = list(map(int, p[1:3])) # start, end
            tss = s + 1 if p[5] == '+' else e
            tssU = tss - tss_up if p[5] == '+' else tss + tss_up
            tssD = tss + tss_down if p[5] == '+' else tss - tss_down
            # fix tss down, within genebody
            tssD = e if tssD > e else s + 1 if tssD < s else tssD
            yield (tssU, tssD, p) if p[5] == '+' else (tssD, tssU, p)


def get_genebody_region(x, gb_up=250, gb_down=2250, tes_extra=0, **kwargs):
    """
    Parameters:
    -----------
    x : str
        gene records in BED6 format
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
    ## is_valid_bed(x)
    with xopen(x) as r:
        for l in r:
            p = l.strip().split('\t') # BED6
            s, e = list(map(int, p[1:3])) # start, end
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


def bed2gtf(x, feature='gene'):
    """
    x : list
        bed record, BED6
    """
    if len(x) < 3:
        return None
    n1 = '{}:{}-{}'.format(x[0], x[1], x[2])
    name, strand = [x[3], x[5]] if len(x) > 5 else [n1, '+']
    s, e = x[1:3] # start, end
    s = int(s) + 1
    des = 'gene_id "{}"; gene_name "{}"'.format(name, name)
    gtf = [x[0], 'bed', feature, s, e, '.', strand, '.', des]
    return list(map(str, gtf))


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
    msg1 = '[TSS-{},TSS+{}]'.format(args['tss_up'], args['tss_down'])
    n1 = args['name']+'.TSSR.{}_TSS_{}'.format(args['tss_up'], args['tss_down'])
    # 1. TSS region
    args['out_dir'] = fix_out_dir(args['out_dir'])
    out = os.path.join(args['out_dir'], n1+'.'+fmt)
    print('- get TSS regions: {}'.format(msg1))
    with xopen(out, 'wt') as w:
        for tss_left, tss_right, p in get_tss_region(x, **args):
            tss = list(map(str, [p[0], tss_left, tss_right]+p[3:6]))
            if fmt == 'gtf':
                tss = bed2gtf(tss, 'gene')
            w.write('\t'.join(tss)+'\n')
    return out


def get_genebody_file(x, fmt='bed', **kwargs):
    """
    x : str
        gene record in BED format, file
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


def count_region(x, region='TSS', **kwargs):
    """
    Count reads on genebody region
    x : str
        gene record in BED format, file
    kwargs : optional arguments
        required, out_dir, gb_up, gb_down, tes_extra, name, overwrite
    """
    args = {
        'bam': None,
        'out_dir': None,
    }
    args.update(kwargs)
    args['out_dir'] = fix_out_dir(args['out_dir'])
    # args['out_dir'] = os.path.join(args['out_dir'],region)  # subdir
    bname = os.path.basename(os.path.splitext(args['bam'])[0]) # bam name
    args['out_dir'] = os.path.join(args['out_dir'], bname)
    if not isinstance(args['name'], str):
        args['name'] = os.path.splitext(os.path.basename(x))[0]
    # 1. extract genebody region
    get_file = get_genebody_file if region == 'genebody' else get_tss_file
    gtf = get_file(x=x, fmt='gtf', **args)
    # 2. count genebody region
    args_fc = {
        'gtf': gtf,
        'bam_list': args['bam'],
        'outdir': args['out_dir'],
        'prefix': os.path.basename(os.path.splitext(gtf)[0])+'.txt',
        'feature_type': 'gene',
    }
    fc = FeatureCounts(**args_fc)
    fc.run()
    return fc.count_txt


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='bed', required=True, 
        help='The gene record in BED6 format')
    parser.add_argument('-b', '--bam', dest='bam', required=True,
        help='bam file, sorted and indexed')
    parser.add_argument('-o', dest='out_dir', required=False, 
        help='The out_dir')
    parser.add_argument('-u', '--tss-up', dest='tss_up', type=int,
        default=150, help='for TSS region, upstream of TSS, default: [150]')
    parser.add_argument('-d', '--tss-down', dest='tss_down', type=int,
        default=150, help='for TSS region, downstream of TSS, default: [150]')
    parser.add_argument('-U', '--gb-up', dest='gb_up', type=int,
        default=250, help='for genebody region, downstream of TSS, default: [250]')
    parser.add_argument('-D', '--gb-down', dest='gb_down',
        default=2250, help='for genebody region, int or "TES", default: [2250]')
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
    cal_pausing_index(**args)
    # tss_f = count_region(x=args.bed, region='TSS', **vars(args))
    # gb_f = count_region(x=args.bed, region='genebody', **vars(args))


if __name__ == '__main__':
    main()

# EOF