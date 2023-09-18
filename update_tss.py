#!/usr/bin/env python
#-*- encoding:utf-8 -*-
"""
Update TSS for gene records (BED)
"""

import os
import argparse


def load_tss(x):
    """
    x is TSS in BED format
    return:
        key = chr:gene 
        value = tss
    """
    d = {}
    with open(x) as r:
        for l in r:
            s = l.strip().split()
            k = f'{s[0]}:{s[3]}' # chr:gene
            d.update({k:s[2]})
    return d


def update_tss(in_bed, tss, out_dir=None, overwrite=False):
    """
    Update the TSS for BED

    Parameters:
    -----------
    in_bed : str
        file, gene records in BED6 format
    tss : str
        file, TSS records in BED6 format, single-base
    out_bed : str, None
        str or None, save the updated gene-records
    """
    print('- Update TSS for genes: {}'.format(os.path.basename(in_bed)))
    # 1. tss
    dt = load_tss(tss) # key=chr:gene, value=tss
    # 2. gene_bed
    if out_dir is None:
        out_bed = os.path.splitext(in_bed)[0]+'.TSS_updated.bed'
    elif isinstance(out_dir, str):
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        out_bed = os.path.join(out_dir, os.path.basename(in_bed))
    else:
        print(f'error: out_dir illegal: {out_dir}')
    # 3. output
    if os.path.exists(out_bed) and overwrite is False:
        print('update_tss() skipped, file exists: {}'.format(out_bed))
        return out_bed
    # 4. run
    with open(in_bed) as r, open(out_bed, 'wt') as w:
        for l in r:
            p = l.strip().split('\t')
            s, e = list(map(int, p[1:3])) # start, end
            tss = dt.get(f'{p[0]}:{p[3]}', None) # key=chr:gene
            if not tss is None:
                tss = int(tss)
                if p[5] == '+':
                    s = tss - 1
                else:
                    e = tss
            if s >= e: # !!! illegal genes, with same id
                print('illegal gene coordinates, {} {}, {}'.format(s, e, p))
                # continue
            else:
                p[1], p[2] = [s, e]
            p = list(map(str, p))
            w.write('\t'.join(p)+'\n')
    print('Saving to file: {}'.format(out_bed))
    return out_bed


# def update_tss(in_bed, tss, out_dir=None, overwrite=False):
#     """
#     Update the TSS for BED

#     Parameters:
#     -----------
#     in_bed : str
#         file, gene records in BED6 format
#     tss : str
#         file, TSS records in BED6 format, single-base
#     out_bed : str, None
#         str or None, save the updated gene-records
#     """
#     print('- Update TSS for genes: {}'.format(os.path.basename(in_bed)))
#     # for genes
#     if out_dir is None:
#         out_bed = os.path.splitext(in_bed)[0]+'.TSS_updated.bed'
#     elif os.path.isdir(out_dir):
#         out_bed = os.path.join(out_dir, os.path.basename(in_bed))
#     if os.path.exists(out_bed) and overwrite is False:
#         print('update_tss() skipped, file exists: {}'.format(out_bed))
#         return out_bed
#     # load TSS
#     t = {}
#     i = 0
#     with open(tss) as r:
#         for l in r:
#             i += 1
#             p = l.split('\t')
#             s,e = list(map(int, p[1:3])) # start, end
#             if e - s > 1:
#                 raise ValueError('line-{}, not a single-base, {}'.format(i, l))
#             t.update({p[3]:p[2]}) # gene_name, TSS
#     # iterate genes
#     with open(in_bed) as r, open(out_bed, 'wt') as w:
#         for l in r:
#             p = l.strip().split('\t')
#             b = p.copy()
#             s, e = list(map(int, b[1:3])) # start, end
#             if s >= e: # !!! illegal genes, with same id !!! 
#                 print('illegal gene coordinates, {} {}, {}'.format(s, e, p))
#                 continue
#             tss = s + 1 if b[5] == '+' else e # current TSS
#             tss_new = t.get(p[3], tss) # updated TSS
#             tss_new = int(tss_new)
#             b[1], b[2] = [tss_new - 1, e] if p[5] == '+' else [s, tss_new]
#             if b[1] >= b[2]: # !!! illegal genes, with same id
#                 print('illegal gene coordinates, {} {}, {}'.format(b[1], b[2], p))
#                 continue
#             b = list(map(str, b))
#             w.write('\t'.join(b)+'\n')
#     # log
#     print('Saving to file: {}'.format(out_bed))
#     return out_bed


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in-bed', dest='in_bed', required=True,
        help='The genebody in BED format,')
    parser.add_argument('-t', '--tss', dest='tss', required=False, 
        help='The TSS record in BED6 format')
    parser.add_argument('-o', '--out-dir', dest='out_dir', required=False,
        default=None, help='The out_dir')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite exists file')
    return parser


def main():
    args = get_args().parse_args()
    update_tss(args.in_bed, args.tss, args.out_dir, args.overwrite)


if __name__ == '__main__':
    main()

# EOF