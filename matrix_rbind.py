#!/usr/bin/env python

"""
Bind multiple matrix, by rows

- merge group names
- keep group names

Usage:
$ python matrix_rbind -m m1.gz m2.gz -o out.gz
"""

import os
import argparse
import json
from xopen import xopen
from utils import load_matrix_header, update_obj, log
from difflib import SequenceMatcher
from functools import reduce

## h1
# '{"upstream": [2000], "downstream": [2000], "body": [2000], "bin size": [100], "ref point": [null],
# "verbose": false, "bin avg type": "mean", "missing data as zero": false, "min threshold": null,
# "max threshold": null, "scale": 1, "skip zeros": false, "nan after end": false, "proc number": 24, 
# "sort regions": "keep", "sort using": "mean", "unscaled 5 prime": [0], "unscaled 3 prime": [0], 
# "group_labels": ["genes"], "group_boundaries": [0, 26799], "sample_labels": ["ChrRNA_mESC_WT.r1_fwd"],
# "sample_boundaries": [0, 60]}'

## h2
# '{"upstream": [2000], "downstream": [2000], "body": [2000], "bin size": [100], "ref point": [null], 
# "verbose": false, "bin avg type": "mean", "missing data as zero": false, "min threshold": null, 
# "max threshold": null, "scale": 1, "skip zeros": false, "nan after end": false, "proc number": 24,
# "sort regions": "keep", "sort using": "mean", "unscaled 5 prime": [0], "unscaled 3 prime": [0],
# "group_labels": ["genes"], "group_boundaries": [0, 26916], "sample_labels": ["ChrRNA_mESC_WT.r1_rev"], 
# "sample_boundaries": [0, 60]}'

# fix groups
class Matrix_rbind(object):
    """
    merge mulitple matrix, by rows
    support only 1 group, n samples
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        
        
    def init_args(self):
        args = {
            'm': None,
            'o': None,
            'overwrite': False
        }
        self = update_obj(self, args, force=False)
        if not isinstance(self.m, list):
            raise ValueError('m, expect list, got {}'.format(
                type(self.m).__name__
            ))
        if not all(list(map(os.path.exists, self.m))):
            raise ValueError('m, file not exists')
        if not isinstance(self.o, str):
            raise ValueError('o, expect str, got {}'.format(
                type(self.o).__name__
            ))
        out_dir = os.path.dirname(os.path.abspath(self.o))
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        
    def fix_group_b(self, x):
        """
        group boundaries
        
        Parameters
        ----------
        x : list of header
        
        merge groups boundaries
        input: [0, 100], [0, 100], ...
        output: [0, 200]
        """
        a = sum([min(i.get('group_boundaries')) for i in x])
        b = sum([max(i.get('group_boundaries')) for i in x])
        return [a, b]


    def fix_group_labels(self, x):
        """
        group labels
        
        Parameters
        ----------
        x : list of header
        
        merge all group labels
        input: ['g1'] ['g1']
        output: ['g1']
        """
        a = list(set(i.get('group_labels')[0] for i in x))
        if len(a) > 1:
            log.error('more than 1 groups found, choose: {}'.format(a[0]))
            a = a[:1]
        return a
    
    
    def fix_sample_labels(self, x):
        """
        sample labels
        
        Parameters
        ----------
        x : list of header
        
        merge samples labels, trim '_fwd', '_rev'
        input: ChrRNA_mESC_WT.r1_fwd, ChrRNA_mESC_WT.r1_rev
        output: ChrRNA_mESC_WT.r1
        """
        g1 = [i.get('sample_labels') for i in x]
        g2 = [[j[i] for j in g1] for i in range(len(g1[0]))]
        s = [reduce(self.find_lcs, i) for i in g2]
        # trim suffix
        return [i.rstrip('._-') for i in s]    
    
    
    def find_lcs(self, s1, s2):
        """
        Find longest common str
        """
        if isinstance(s1, str) and isinstance(s2, str):
            m = SequenceMatcher(None, s1, s2) # match
            l = m.find_longest_match(0, len(s1), 0, len(s2))
            out = s1[l.a:(l.a+l.size)]
        else:
            out = None
        return out
    
    
    def match_val(self, s1, s2):
        """
        a and b are identical
        """
        return s1 == s2
     
    
    def run(self):
        if os.path.exists(self.o) and not self.overwrite:
            log.info('matrix_rbind() skipped, file exists: {}'.format(self.o))
        else:
            h = [load_matrix_header(i) for i in self.m]
            h1 = h[0] # first one
            args = {
                'group_labels': self.fix_group_labels(h),
                'group_boundaries': self.fix_group_b(h),
                'sample_labels': self.fix_sample_labels(h),
            }
            q = [[p.get(i) for p in h] for i,j in h1.items() if i not in args]
            t = [reduce(self.match_val, i) for i in q]
            if not all(t):
                raise ValueError('parameters not matched')
            # update header
            h1.update(args)
            hs = json.dumps(h1) #
            # hs = '@'+hs # add @
            # merge all lines
            with xopen(self.o, 'wt') as w:
                w.write('@'+hs+'\n') # header, new
                for hx in self.m:
                    with xopen(hx) as r:
                        for line in r:
                            if line.startswith('@'):
                                continue
                            w.write(line)     


def get_args():
    example = ' '.join([
        '$ plotProfile',
        '-m input.mat -o out.png',
        '--samplesLabel A B C --regionsLabel gene1 gene2',
        '--colors black lightblue --yMin 0 --yMax 0.4',
        '--perGroup',
    ])
    parser = argparse.ArgumentParser(prog='bw2matrix.py',
                                     description='bw2matrix',
                                     epilog=example,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-m', dest='m', nargs='+', required=True,
        help='matrix files, by computeMatrix ')
    parser.add_argument('-o', dest='o', required=True,
        help='output matrix')
    parser.add_argument('--overwrite', dest='overwrite', action='store_true',
        help='Overwrite output file')
    return parser


def main():
    args = vars(get_args().parse_args())
    Matrix_rbind(**args).run()

    
if __name__ == '__main__':
    main()

# EOF