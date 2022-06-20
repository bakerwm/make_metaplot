#!/usr/bin/env python
#-*- encoding: utf-8 -*-
"""
Drecated, see `cal_pausing_index.py` 

Pick the strongest TSS for gene, with multiple TSSs
criteria
- the strongest RNAP2 ChIP-seq signal (TSS region: tss_up + tss_down)
- the distal TSS, if multiple TSShave equal ChIP-seq sigals
- H3K4me3 enrichment (ChIP/input) > 4 (optional) 

input files:
1. TSS annotation, (from ENSEMBL-BioMart, UCSC-TableBrowser,...), BED6 
2. RNAP2 ChIP-seq (BAM), (Rbp1, ..., NOT S5P, S2P, ...)


output files
1. TSS region (BED6), single-base
"""


import os
import sys
import re
import pathlib
import argparse
from xopen import xopen
from hiseq.utils.featurecounts import FeatureCounts
from hiseq.utils.bam import Bam
# from hiseq.utils.utils import Config
# from logging import raiseExceptions


def pick_tss(x, bam, **kwargs):
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
    print('> Pick the strongest TSS')
    args = {
        'out_dir': None,
        'tss_up': 150,
        'tss_down': 150,
        'name': None,
        'overwrite': False,
    }
    args.update(kwargs)
    # 1. prepare files
    args['out_dir'] = fix_out_dir(args['out_dir'])
    if not isinstance(args['name'], str):
        args['name'] = os.path.basename(os.path.splitext(x)[0])
    prefix = os.path.join(args['out_dir'], args['name']) 
    tss_r = prefix+'.TSSR.{}_TSS_{}.bed'.format(args['tss_up'], args['tss_down'])
    # 2. extract TSS region
    get_tss_region(x, tss_r, args['tss_up'], args['tss_down'])
    # 3. count reads on TSS regions
    args['name'] = os.path.basename(os.path.splitext(tss_r)[0])
    tss_fc = count_region(x=tss_r, bam=bam, **args)
    # 4. pick TSS sites
    tss_fix = os.path.join(os.path.dirname(tss_fc), args['name']+'.TSS_fixed.bed')
    pick_distal_tss(tss_fc, tss_fix, args['overwrite'])
    # 5. convert TSS to single-base BED format
    tss_fix_new = recover_tss(tss_fix) # or retrieve from total TSS, by gene_name
    return tss_fix_new


def pick_distal_tss(x, y, overwrite=False):
    """
    Pick the distal TSS for genes with multiple TSSs that equal in RNAP2 signal
    choose if RNAP2 signal was zero

    x : str
        otuput of featureCounts, read counts in column-7
    """
    if os.path.exists(y) and overwrite is False:
        print('pick_distal_tss() skipped, file exists.')
        return None
    d = {}
    for p in load_fc(x):
        # if float(p[-1]) == 0:
        #     continue
        g,_ = p[0].split(':', 1) # fix gene_name:tss
        dg = d.get(g, [])
        dg.append(p)
        d.update({g:dg})
    # pick strongest tss
    print('Number of TSSs for gene: {}'.format(len(d)))
    # choose distal tss, for multiple TSSs with equal counts
    with open(y, 'wt') as w:
        for g,fc in d.items():
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


def recover_tss(x):
    """
    Recover the TSS region to TSS site
    parse tss_up and tss_down from filename: ...TSSR.150_TSS_150.TSS_fixed.bed
    
    Parameters:
    -----------
    x : str
        file, fixed TSS in BED format, see pick_distal_tss()
        suffix: ".TSSR.150_TSS_150.TSS_fixed.bed"
    
    write the final results to file: ".TSS_fixed.bed"
    """
    if not os.path.exists(x):
        print('file not exists: {}'.format(x))
        sys.exit(1)
    # check file name
    p = re.compile('(\d+)_TSS_(\d+).TSS_fixed')
    g = p.search(x)
    if g is None:
        print('x, expect [.150_TSS_150.TSS_fixed], got {}'.format(x))
        sys.exit(1)
    tss_up, tss_down = (int(g.group(1)), int(g.group(2)))
    # new name
    x_new = re.sub('TSSR.\d+_TSS_\d+\.', '', x, flags=re.IGNORECASE)
    i = 0
    with open(x_new, 'wt') as w, open(x) as r:
        for l in r:
            i += 1
            p = l.strip().split('\t')
            s,e = list(map(int, p[1:3])) # start, end
            p[1], p[2] = (s+tss_up, e-tss_down) if p[5] == '+' else (s+tss_down, e-tss_up)
            if p[1] >= p[2]:
                p[1] = p[2] - 1 # one the edges
                print('line-{}, illegal BED, {}'.format(i, x))
            b = list(map(str, p))
            w.write('\t'.join(b)+'\n')
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


def get_tss_region(x, y, tss_up=0, tss_down=0):
    """
    Parameters:
    -----------
    x : str
        TSS records in BED format, tss_start, tss_end, single base
    y : str
        output file, formated TSS record in BED6 
    tss_up : int
        TSS region, upstream of TSS, default: [0]
    tss_down : int
        TSS region, downstream of TSS, default: [0]
    # format TSS 
    1. slop TSS region: tss_up, tss_down
    # return
    chr,tss-1,tss,gene_name,254,strand
    # depict, how to define the TSS region
    # genes on forward strand(+):
                  TSS                                        TES    
    (+)------------|----------------[gene]--------------------|--------->
        |--tss_up--|--tss_down--|
        (.......TSS rgion.......)

    ## genes on reverse strand(-):
                  TES                                        TSS    
    (-)<-----------|----------------[gene]--------------------|-----------
                                                 |--tss_down--|--tss_up--|
                                                 (.......TSS rgion.......)
    """
    with xopen(x) as r, xopen(y, 'wt') as w:
        for l in r:
            p = l.strip().split('\t')
            s,e = list(map(int, p[1:3]))
            # update start,end
            s,e = (s-tss_up, e+tss_down) if p[5] == '+' else (s-tss_down, e+tss_up)
            # fix start
            # fix chr size, ! to-do
            if s < 0:
                s = 0
            b = [p[0], s, e] + p[3:6]
            b = list(map(str, b))
            if y.endswith('.gtf'):
                b = bed2gtf(b, feature='gene')
            w.write('\t'.join(b)+'\n')


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


## count reads on TSS region
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


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='bed', required=True, 
        help='The TSS record in BED6 format, see "format_tss.py" for help')
    parser.add_argument('-b', '--bam', dest='bam', required=True,
        help='bam file, sorted and indexed')
    parser.add_argument('-o', dest='out_dir', required=False, 
        help='The out_dir')
    parser.add_argument('-u', '--tss-up', dest='tss_up', type=int,
        default=150, help='for TSS region, upstream of TSS, default: [150]')
    parser.add_argument('-d', '--tss-down', dest='tss_down', type=int,
        default=150, help='for TSS region, downstream of TSS, default: [150]')
    parser.add_argument('-n', '--name', dest='name', type=str,
        default=None, help='name of the output files, default: [auto]')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite exists file')
    return parser


def main():
    args = vars(get_args().parse_args())
    bed = args.pop('bed')
    bam = args.pop('bam')
    pick_tss(bed, bam, **args)


if __name__ == '__main__':
    main()

# EOF