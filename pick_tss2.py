#!/usr/bin/env python 
"""
Pick the TSS position, identified by multiple experiments

input: 
  - list of BED file, TSS (single base) 

output:
  - BED file
"""

import os
import sys
import argparse


def load_bed(x):
    """
    load TSS BED file, name:bed as dict
    """
    d = {}
    with open(x) as r:
        for l in r:
            p = l.strip().split('\t')
            id = '{}:{}'.format(p[0], p[3]) # chr:gene
            d.update({id:p})
    return d


def load_bed2(x):
    """
    Parameters:
        x : list, BED files
    load list of TSS BED files, dicts
    """
    return [load_bed(i) for i in x]


def most_freq(x):
    """
    choose the most frequent value in x (list)
    """
    return max(set(x), key = x.count)


def pick_tss2(x, out):
    """
    Parameters:
        x : list, BED files, single base    
    """
    dn = load_bed2(x)
    d1 = dn.pop() # first element
    with open(out, 'wt') as w:
        for k,v in d1.items():
            t0 = [i.get(k, None) for i in dn]
            t1 = [i[1] for i in t0 if not i is None]
            t1.append(v[1]) #
            v[1] = most_freq(t1) # update TSS-1
            v[2] = str(int(v[1]) + 1)
            w.write('\t'.join(v)+'\n')


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tss', required=True, nargs='+',
        help='TSS for each experiment')
    parser.add_argument('-o', dest='out_bed', required=True, 
        help='The out_bed')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite exists file')
    return parser


def main():
    args = get_args().parse_args()
    # tss names 
    bam_names = [os.path.basename(os.path.dirname(i)) for i in args.tss]
    msg = 'selected {} experiments: [{}]'.format(len(bam_names), ','.join(bam_names))
    print(msg)
    pick_tss2(args.tss, args.out_bed)


if __name__ == '__main__':
    main()

# EOF