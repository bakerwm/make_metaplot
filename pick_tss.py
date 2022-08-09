#!/usr/bin/env python
#-*- encoding: utf-8 -*-
"""
Update TSS for genes by Pol II ChIP data

criteria: pick the "TSS" with the strongest Pol II ChIP signal
arguments: gene_bed, tss_bed, bam, ...
"""


from operator import truediv
import os
import sys
import re
import pathlib
import argparse
import pysam
from datetime import datetime
from xopen import xopen
from hiseq.utils.featurecounts import FeatureCounts
# from hiseq.utils.file import symlink_file


def pick_tss(x, bam, **kwargs):
    """
    Pick Top/Distal TSS

    Parameters:
    -----------
    x : str
        TSS records in BED6 format, see script "format_tss.py" for help
    bam : str
        The RNAP2 ChIP alignment file, to determine TSS site
    out_dir : str
        Directory to save final records
    tss_up : int
        TSS region, upstream of TSS, default: [150]
    tss_down : int
        TSS region, downstream of TSS, default: [150]
    kwargs
    """
    args = {
        'prefix': None,
        'overwrite': False
    }
    args.update(kwargs)
    overwrite = args.get('overwrite', False)
    # 1. extract TSS region
    tss_f = get_tss_file(x, fmt='bed', **args)
    # 2. count reads on TSS regions
    args['prefix'] = os.path.basename(os.path.splitext(tss_f)[0]) # update name
    tss_fc = count_region(x=tss_f, bam=bam, **args)
    # 3. pick TSS sites
    tss_top = pick_top_tss(tss_fc, overwrite)
    # 4. convert TSS to single-base BED format
    tss_top_fixed = recover_tss(tss_top) # retrieve from total TSS, by gene_name
    # # 5. save top_tss to main directory
    # tss = os.path.join(args['out_dir'], os.path.basename(tss_top_fixed))
    # symlink_file(tss_top_fixed, tss)
    return tss_top_fixed #_top_fixed


def pick_top_tss(x, overwrite=False):
    """
    Pick Top/Distal TSS

    Parameters:
    -----------
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
    # choose top tss
    # choose distal tss, for multiple TSSs with equal counts
    with open(x_new, 'wt') as w:
        for k,fc in d.items(): # k=gene+chr
            g,c = k.split(':', 1)
            n = [float(i[6]) for i in fc] # count in column-7
            max_idx = [i for i,j in enumerate(n) if j == max(n)] # top TSS
            if len(max_idx) == 1:
                max_fc = fc[max_idx[0]]
            else:
                fc2 = [j for i,j in enumerate(fc) if i in max_idx]
                max_fc = pick_distal_tss(fc2) # distal TSS
            # fix BED start, 0-index
            max_fc[2] = str(int(max_fc[2])-1) 
            t = max_fc[1:4]+[g, 254, max_fc[4]] # BED format
            t = list(map(str, t))
            w.write('\t'.join(t)+'\n')
    return x_new


def pick_distal_tss(x):
    """
    Parameters:
    -----------
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
    g = re.search('\.TSSR\.(\d+)_TSS_(\d+)', x)
    if g is None:
        print('x, expect [.TSSR.150_TSS_150], got {}'.format(x))
        sys.exit(1)
    tss_up, tss_down = (int(g.group(1)), int(g.group(2)))
    # new name
    x_new = re.sub('.TSSR.\d+_TSS_\d+\.top_TSS', '.TSS', x, flags=re.IGNORECASE)
    if os.path.exists(x_new) and overwrite is False:
        print('recover_tss() skipped, file exists: {}'.format(x_new))
        return x_new
    i = 0
    with open(x_new, 'wt') as w, open(x) as r:
        for l in r:
            i += 1
            p = l.strip().split('\t')
            s,e = list(map(int, p[1:3])) # start, end
            s,e = (s+tss_up, e-tss_down) if p[5] == '+' else (s+tss_down, e-tss_up)
            if s >= e:
                s = e - 1 # on the edges
                print('line-{}, illegal BED, {}, {}'.format(i, s, e))
            p[1], p[2] = list(map(str, [s, e]))
            w.write('\t'.join(p)+'\n')
    return x_new


def get_tss_file(x, fmt='bed', **kwargs):
    """
    Parameters:
    -----------
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
        'prefix': None,
        'overwrite': False,
    }
    args.update(kwargs)
    # update files
    if not isinstance(args['prefix'], str):
        args['prefix'] = os.path.basename(os.path.splitext(x)[0])
    args['out_dir'] = fix_out_dir(args['out_dir'])
    n = args['prefix']+'.TSSR.{}_TSS_{}'.format(args['tss_up'], args['tss_down'])
    out = os.path.join(args['out_dir'], n+'.'+fmt)
    # message
    msg = 'TSS region: --[({}bp)--TSS--({}bp)]--'.format(args['tss_up'], args['tss_down'])
    print(msg)
    if os.path.exists(out) and args.get('overwrite', False) is False:
        print('file exists: {}'.format(out))
    else:
        with xopen(out, 'wt') as w:
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
    # print('!B-2', kwargs, type(tss_up), type(tss_down))
    # is_valid_bed(x)
    with xopen(x) as r:
        for l in r:
            p = l.strip().split('\t') # BED6
            if len(p) < 6:
                continue
            s,e = list(map(int, p[1:3])) # start, end
            s,e = (s-tss_up, s+tss_down+1) if p[5] == '+' else (e-tss_down-1, e+tss_up)
            if s < 0:
                s = 0 # fix border
            b = [p[0], s, e] + p[3:6]
            yield list(map(str, b)) # return BED record


# count reads on TSS region
def count_region(x, bam, **kwargs):
    """
    Count reads on regions using featureCounts

    Parameters:
    -----------
    x : str
        gene record in BED/GTF format, file
    bam : str
        Alignment file of Pol II ChIP-seq
    kwargs : optional arguments
        required, out_dir, name, overwrite
    """
    args = {
        'out_dir': None,
        'prefix': None,
        'overwrite': False,
    }
    args.update(kwargs)
    # update files
    args['out_dir'] = fix_out_dir(args['out_dir'])
    bname = os.path.basename(os.path.splitext(bam)[0]) # bam name
    args['out_dir'] = os.path.join(args['out_dir'], bname) # update subdir
    if not isinstance(args['prefix'], str):
        args['prefix'] = os.path.splitext(os.path.basename(x))[0]
    # 1. gtf file
    if x.endswith('.gtf'):
        gtf = x
    elif x.endswith('.bed'):
        gtf = bed2gtf_file(x)
    else:
        sys.exit('unknown x={}'.format(x))
    # 2. count genebody region
    args_fc = {
        'gtf': gtf,
        'bam_list': bam,
        'outdir': args['out_dir'],
        'prefix': args['prefix']+'.txt',
        'feature_type': 'gene',
    }
    # index bam file
    if not os.path.exists(bam+'.bai'):
        pysam.index(bam)
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


def bed2gtf_file(x, y=None, feature='gene', overwrite=False):
    """
    Convert BED file to GTF file
    """
    # output file
    if isinstance(y, str):
        y = os.path.abspath(y)
        out_dir = os.path.dirname(y)
    else:
        x = os.path.abspath(x)
        out_dir = os.path.dirname(x)
        y = os.path.splitext(x)[0] + '.gtf'
    out_dir = fix_out_dir(out_dir)
    # run
    if os.path.exists(y) and overwrite is False:
        print('file exists: {}'.format(y))
    else:
        with open(x) as r, open(y, 'wt') as w:
            for l in r:
                p = l.strip().split('\t')
                g = bed2gtf(p, feature='gene')
                w.write('\t'.join(g)+'\n')
    return y


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tss-bed', dest='tss_bed', required=True, 
        help='The TSS record in BED6 format, see "format_tss.py" for help')
    parser.add_argument('-b', '--tss-bam', dest='tss_bam', required=False,
        help='The alignemnt file of Pol2 ChIP, in BAM format')
    parser.add_argument('-o', dest='out_dir', required=False, 
        help='The out_dir')
    parser.add_argument('-u', '--tss-up', dest='tss_up', type=int,
        default=150, help='for TSS region, upstream of TSS, default: [150]')
    parser.add_argument('-d', '--tss-down', dest='tss_down', type=int,
        default=150, help='for TSS region, downstream of TSS, default: [150]')
    parser.add_argument('-n', '--name', dest='prefix', type=str,
        default=None, help='name of the output files, default: [auto]')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite exists file')
    return parser


def main():
    args = vars(get_args().parse_args())
    # gene = args.pop('gene_bed')
    tss_bed = args.pop('tss_bed', None)
    tss_bam = args.pop('tss_bam', None)
    pick_tss(tss_bed, tss_bam, **args)


if __name__ == '__main__':
    main()

# EOF