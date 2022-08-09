#!/usr/bin/env python
#-*- encoding: utf-8 -*-
"""
Pick ChIP occupied genes

ChIP peaks overlapped on TSS of genes
input: 
1. gene - bed file
2. peak - bed file
3. type - tss, genebody
output:
1. gene - bed file
"""

import os
import sys
import pathlib
import argparse
import pybedtools


def pick_tss(x, up=0, down=0):
    with open(x) as r:
        for l in r:
            s = l.strip().split('\t')
            if len(s) < 6: 
                continue # BED6
            s[1], s[2] = [s[1], int(s[1])+1] if s[5] == '+' else [int(s[2])-1, s[2]]
            # update flanking
            s[1], s[2] = [int(s[1])-up, int(s[2])+down] if s[5] == '+' else [int(s[1])-down, int(s[2])+up]
            yield list(map(str, s[:6]))


def get_tss_file(x, out=None, up=0, down=0):
    if not isinstance(out, str):
        out = os.path.splitext(x)[0]+'.tss.bed'
    if os.path.exists(out):
        print('file exists: {}'.format(out))
        return out
    with open(out, 'wt') as w:
        for b in pick_tss(x, up, down):
            w.write('\t'.join(b)+'\n')
    return out


def pick_chip_gene(gene, peak, out=None, on_tss=False, up=0, down=0):
    """
    genes overlapped with peaks
    bedtools intersect -u -a ${b} -b ${pol2_peak} > ${out}
    """
    if not isinstance(out, str):
        out = os.path.splitext(gene)[0] + '.peak.bed'
    out = os.path.abspath(out)
    out_dir = os.path.dirname(out)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # 1. on genebody
    if not on_tss:
        pybedtools.BedTool(gene).intersect(u=True, b=peak).saveas(out)
    else:
        # 2. on TSS
        ### pick tss bed
        tss = os.path.splitext(out)[0] + '.tss.bed'
        get_tss_file(gene, tss, up, down)
        tb = pybedtools.BedTool(tss).intersect(u=True, b=peak)
        tg = {d.name:1 for d in tb}
        ### recover to gene
        with open(gene) as r, open(out, 'wt') as w:
            for l in r:
                s = l.strip().split('\t')
                if s[3] in tg:
                    w.write('\t'.join(s)+'\n')
        ### remove tss file
        if os.path.exists(tss):
            os.remove(tss)
    # output
    print('save to file: {}'.format(out))


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gene-bed', dest='gene_bed', required=True,
        help='The genebody in BED format,')
    parser.add_argument('-p', '--peak', dest='peak_bed', required=True, 
        help='The peak record in BED format')
    parser.add_argument('-o', dest='out', required=False, default=None, 
        help='The output file')
    parser.add_argument('-u', '--up', dest='up', type=int,
        default=0, help='for TSS region, upstream of TSS, default: [0]')
    parser.add_argument('-d', '--down', dest='down', type=int,
        default=0, help='for TSS region, downstream of TSS, default: [0]')
    parser.add_argument('--tss', dest='on_tss', action='store_true',
        help='Peaks overlapped on TSS of genes')
    return parser


def main():
    args = get_args().parse_args()
    pick_chip_gene(args.gene_bed, args.peak_bed, args.out, args.on_tss, args.up, args.down)


if __name__ == '__main__':
    main()

# EOF