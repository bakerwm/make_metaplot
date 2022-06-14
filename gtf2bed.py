#!/usr/bin/env python 
#-*- encoding: UTF-8 -*-
"""
Extract gene from GTF of ENSEMBL
"""

import os
import sys
import io
import pathlib
import argparse
import pybedtools
from datetime import datetime
from xopen import xopen


def gtf_to_bed(gtf, bed, **kwargs):
    """
    Parameters:
    -----------
    feature : str
        gene|exon|CDS|transcript, default: None
    gene_biotype : str
        TEC|protein_coding, default: None
    Convert GTF to BED format
    filt by: gene_biotype (ENSEMBL), gene_type (GENCODE)
    """
    overwrite = kwargs.get('overwrite', False)
    if os.path.exists(bed) and overwrite is False:
        print('file exists, gtf_to_bed() skipped, {}'.format(bed))
    else:
        with xopen(gtf) as r, xopen(bed, 'wt') as w:
            for g in read_gtf(r, **kwargs):
                g_chr = g[0] if g[0].startswith('chr') else 'chr'+g[0]
                g_start = str(int(g[3])-1)
                g_name = parse_gtf_desc(g[8], 'gene_name')
                if g_name is None:
                    g_name = parse_gtf_desc(g[8], 'gene_id') # alternative
                gene_biotype = parse_gtf_desc(g[8], 'gene_biotype')
                bed6 = [g_chr, g_start, g[4], g_name, '254', g[6]]
                if kwargs.get('append_biotype', False):
                    bed6.append(gene_biotype)
                bed6 = list(map(str, bed6))
                w.write('\t'.join(bed6)+'\n')


def read_gtf(fh, **kwargs):
    """
    Parameters:
    -----------
    feature : str
        gene|exon|CDS|transcript, default: None (all)
    gene_biotype : str
        TEC|protein_coding, default: None (all)

    filt GTF records by [feature] and [gene_biotype]
    """
    feature = kwargs.get('feature', None)
    gene_biotype = kwargs.get('gene_biotype', None)
    for l in fh:
        p = l.strip().split('\t')
        # filtering
        if len(p) < 9:
            continue
        # filt by feature
        if not feature is None:
            if isinstance(feature, str) and not p[2] == feature:
                continue
            elif isinstance(feature, list) and not p[2] in feature:
                continue
            else:
                pass
        # filt by gene_biotype
        gb = parse_gtf_desc(p[8], 'gene_biotype')
        if not gene_biotype is None:
            if isinstance(gene_biotype, str) and not gb == gene_biotype:
                continue
            elif isinstance(gene_biotype, list) and not gb in gene_biotype:
                continue
            else:
                pass
        yield p


def parse_gtf_desc(s, key='gene_name'):
    """
    Description in GTF:
    gene:
    ENSEMBL:
    gene_id "ENSMUSG00000102693"; gene_version "1"; gene_name "4933401J01Rik"; 
    gene_source "havana"; gene_biotype "TEC"; havana_gene "OTTMUSG00000049935";
    havana_gene_version "1";
    GENCODE:
    gene_id "ENSG00000227232.4"; transcript_id "ENSG00000227232.4"; 
    gene_type "pseudogene"; gene_status "KNOWN"; gene_name "WASH7P"; 
    transcript_type "pseudogene"; transcript_status "KNOWN"; 
    transcript_name "WASH7P"; level 2; havana_gene "OTTHUMG00000000958.1";
    """
    d = {}
    if isinstance(s, str):
        for p in s.strip().split(';'):
            if " " in p:
                a, b = p.split()[:2]
                b = b.replace('"', '')
                d.update({a:b})
    # fix gene_type (GENCODE), gene_biotype (ENSEMBL)
    if key == 'gene_type' and 'gene_biotype' in s:
        key = 'gene_biotype'
    elif key == 'gene_biotype' and 'gene_type' in s:
        key = 'gene_type'
    return d.get(key, None)


def filt_bed_size(x, y, size=5000, overwrite=False):
    """
    gene_biotype saved in column-7 (last column)
    """
    if os.path.exists(y) and overwrite is False:
        print('file exists, flit_bed_size() skipped, {}'.format(y))
    else:
        bed = pybedtools.BedTool(x).filter(lambda i: len(i) >= size) # gene size
        bed.saveas(y)


def filt_by_flank(x, y, flank_size=6000, overwrite=False):
    """
    Filt the BED file by flanking size
    1. flanking BED (6kb) not including genes
    """
    if os.path.exists(y) and overwrite is False:
        print('file exists, flit_by_flank() skipped, {}'.format(y))
        return None
    # sort bed
    xb = pybedtools.BedTool(x).sort()
    # filter by flank size
    with xopen(y, 'wt') as w:
        pre_gene = [None, None] #fwd, rev
        for i in xb:
            if i.strand == '+':
                if pre_gene[0] is None:
                    pre_gene[0] = i
                else:
                    if pre_gene[0].chrom == i.chrom and i.start - pre_gene[0].end <= flank_size:
                        continue
                    else:
                        pre_gene[0] = i
            elif i.strand == '-':
                if pre_gene[1] is None:
                    pre_gene[1] = i
                else:
                    if pre_gene[1].chrom == i.chrom and i.start - pre_gene[1].end <= flank_size:
                        continue
                    else:
                        pre_gene[1] = i
            bed6 = list(map(str, i[:7])) # bed6+1
            w.write('\t'.join(bed6)+'\n')


def filt_gene_biotype(x, y, gene_biotype=None, overwrite=False):
    """
    gene_biotype saved in column-7 (last column)
    """
    if os.path.exists(y) and overwrite is False:
        print('file exists, flit_bed_size() skipped, {}'.format(y))
    else:
        bed = pybedtools.BedTool(x)
        if isinstance(gene_biotype, str):
            bed = bed.filter(lambda i: i[-1] == gene_biotype)
        elif isinstance(gene_biotype, list):
            bed = bed.filter(lambda i: i[-1] in gene_biotype)
        else:
            pass
        bed.saveas(y)


def cur_time():
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='gtf', required=True, 
        help='The gene GTF file')
    parser.add_argument('-o', dest='out_dir', required=False, 
        help='The out_dir')
    parser.add_argument('-f', '--feature', 
        help='[gene|exon|CDS|transcript], default [None] all')
    parser.add_argument('-g', '--gene-biotype', dest='gene_biotype',
        help='The gene_biotype, [protein_coding|...], default: [None] all')
    parser.add_argument('-x', '--append-biotype', dest='append_biotype', 
        action='store_true', help='Add gene_biotype to column-7 in BED')
    parser.add_argument('-gs', '--gene-size', dest='gene_size', type=int, 
        default=500, help='The minimum size of BED feature, default [0]')
    parser.add_argument('-fs', '--flank-size', dest='flank_size', type=int, 
        default=1, help='The minimum distance to neighbour gene, default [0]')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite exists file')
    return parser


def main():
    """
    1. extract all gene
    2. genebody size (750)
    3. flanking size 
    4. protein coding
    """
    args = get_args().parse_args()
    gene_size = args.gene_size
    flank_size = args.flank_size
    out_dir = args.out_dir
    # out_dir
    if not isinstance(out_dir, str):
        out_dir = pathlib.Path.cwd()
    out_dir = os.path.abspath(out_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # prepare files
    gs = '{:d}k'.format(int(gene_size/1000)) if gene_size > 1000 else gene_size
    fs = '{:d}k'.format(int(flank_size/1000)) if flank_size > 1000 else flank_size
    prefix = os.path.splitext(os.path.basename(args.gtf))[0]
    prefix = os.path.splitext(prefix)[0] if prefix.endswith('.gtf') else prefix
    bed_gene = os.path.join(out_dir, prefix+'.gene.bed')
    bed_size = os.path.join(out_dir, prefix+'.gene.g{}.bed'.format(gs))
    bed_pc = os.path.join(out_dir, prefix+'.gene.g{}.pc.bed'.format(gs))
    bed_flank = os.path.join(out_dir, prefix+'.gene.g{}.pc.f{}.bed'.format(gs,fs))

    # step1. gtf to bed, feature, (all genes)
    print('[{}] 1. extract all gene'.format(cur_time()))
    kwargs = {
        'feature': args.feature, 
        'gene_biotype': None, # args.gene_biotype,
        'append_biotype': True, # args.append_biotype,
        'overwrite': args.overwrite,
        }
    gtf_to_bed(args.gtf, bed_gene, **kwargs) # protein_coding, genes

    # step2. filt by gene size
    print('[{}] 2. filter by, gene_size={}'.format(cur_time(), gene_size))
    filt_bed_size(bed_gene, bed_size, size=args.gene_size, overwrite=args.overwrite)

    # step3. filt by genetype (protein_coding)
    print('[{}] 3. filter by, gene_biotype={}'.format(cur_time(), args.gene_biotype))
    filt_gene_biotype(bed_size, bed_pc, gene_biotype=args.gene_biotype, overwrite=args.overwrite)

    # step4. filt by flanking size
    print('[{}] 4. filter by, flank_size={}'.format(cur_time(), flank_size))
    filt_by_flank(bed_pc, bed_flank, flank_size=flank_size, overwrite=args.overwrite)

    print('Finish, save to file: {}'.format(bed_pc))

if __name__ == '__main__':
    main()
    
# EOF