#!/usr/bin/env python
#-*- encoding:utf-8 -*-
"""
Update TSS for gene records (BED)
"""

import os
import argparse


def update_tss(in_bed, in_tss, out_dir=None, overwrite=False):
    """
    Update the TSS for BED

    Parameters:
    -----------
    in_bed : str
        file, gene records in BED6 format
    in_tss : str
        file, TSS records in BED6 format, single-base
    out_bed : str, None
        str or None, save the updated gene-records
    """
    print('- Update TSS for genes: {}'.format(os.path.basename(in_bed)))
    # for genes
    if out_dir is None:
        out_bed = os.path.splitext(in_bed)[0]+'.TSS_updated.bed'
    elif os.path.isdir(out_dir):
        out_bed = os.path.join(out_dir, os.path.basename(in_bed))
    if os.path.exists(out_bed) and overwrite is False:
        print('update_tss() skipped, file exists: {}'.format(out_bed))
        return out_bed
    # load TSS
    t = {}
    i = 0
    with open(in_tss) as r:
        for l in r:
            i += 1
            p = l.split('\t')
            s,e = list(map(int, p[1:3])) # start, end
            if e - s > 1:
                raise ValueError('line-{}, not a single-base, {}'.format(i, l))
            t.update({p[3]:p[2]}) # gene_name, TSS
    # iterate genes
    with open(in_bed) as r, open(out_bed, 'wt') as w:
        for l in r:
            p = l.strip().split('\t')
            b = p.copy()
            s, e = list(map(int, b[1:3])) # start, end
            if s >= e: # !!! illegal genes, with same id !!! 
                print('illegal gene coordinates, {} {}, {}'.format(s, e, p))
                continue
            tss = s + 1 if b[5] == '+' else e # current TSS
            tss_new = t.get(p[3], tss) # updated TSS
            tss_new = int(tss_new)
            b[1], b[2] = [tss_new - 1, e] if p[5] == '+' else [s, tss_new]
            if b[1] >= b[2]: # !!! illegal genes, with same id
                print('illegal gene coordinates, {} {}, {}'.format(b[1], b[2], p))
                continue
            b = list(map(str, b))
            w.write('\t'.join(b)+'\n')
    # log
    print('Saving to file: {}'.format(out_bed))
    return out_bed


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in-bed', dest='in_bed', required=True,
        help='The genebody in BED format,')
    parser.add_argument('-t', '--in-tss', dest='in_tss', required=False, 
        help='The TSS record in BED6 format, see "format_tss.py" for help')
    parser.add_argument('-o', dest='out_dir', required=False, default=None,
        help='The out_dir')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite exists file')
    return parser


def main():
    args = get_args().parse_args()
    if isinstance(args.out_dir, str):
        try:
            if not os.path.exists(args.out_dir):
                os.makedirs(args.out_dir)
        except:
            print('failed, could not create dir: {}'.format(args.out_dir))
            args.out_dir = None
    update_tss(args.in_bed, args.in_tss, args.out_dir, args.overwrite)


if __name__ == '__main__':
    main()

# EOF