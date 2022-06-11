#!/usr/bin/env python
"""
!!! ATTENTION !!!
This script NOT working for general `rbind` matrix.
This script was designed to merge two matrix of sens + anti stranded.
sample_labels contain: '_fwd', '_rev' suffix,
group_labels contain: '_fwd.bed', '_rev.bed' suffix

Requirements:
1. same sample labels
2. same group labels

Why?
merge sense and antisense matrix, into one single matrix.

How?
1. update labels (sample/group)
2. bind rows
3. update group boundaries (for group)


Bind multiple matrix, by rows
row: groups (regions)
column: samples (signal)

update:
- sample_labels
- group_labels
- group_boundaries

Usage:
$ python matrix_rbind -m m1.gz m2.gz -o out.gz
"""

from multiprocessing.sharedctypes import Value
import os
import argparse
import json
from xopen import xopen
from utils import load_matrix, update_obj, log
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
    row: groups (regions)
    column: samples (signal)
    merge mulitple matrix, by rows
    support only 1 group, n samples
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        self.headers = self.load_matrix2(self.m, header_only=True)


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


    def load_matrix2(self, x, header_only=True):
        """
        load header/body from multiple matrix files
        Parameters:
        -----------
        x : list
            list of matrix files
        header_only: bool
            load header only
        """
        y = [load_matrix(i, header_only) for i in x]
        if any([len(i) == 0 for i in y]):
            s = ['{} : {}'.format(i, j) for i,j in zip([len(i) > 0 for i in y], self.m)]
            print('\n'.join(s))
            raise ValueError('illegal matrix file(s) detected')
        return y # list of dict


    def check_labels(self):
        """
        Check sample/group labels
        """
        # 1. sample_labels
        s1 = [i.get('sample_labels', []) for i in self.headers]
        s1x = self.update_labels(s1)
        k1 = self.is_valid_labels(s1x)
        # 2. group labels
        s2 = [i.get('group_labels', []) for i in self.headers]
        s2x = self.update_labels(s2)
        k2 = self.is_valid_labels(s2x)
        # output
        if not all([k1, k2]):
            raise ValueError('group/sample labels not valid')
        return (s1x[0], s2x[0])


    def is_valid_labels(self, x):
        """
        check lables
        - identical in all (length, content)
        """
        k1 = len(set([len(i) for i in x])) == 1
        x2 = [[j[i] for j in x] for i,_ in enumerate(x[0])]
        k2 = all([len(set(i))==1 for i in x2])
        return k1 and k2


    def update_labels(self, x):
        """
        update sample/group lables remove ('_fwd', '_rev')
        input: [wt_fwd, mut_fwd], [wt_rev, mut_rev]
        output: [wt, mut]
        """
        # s = [i.get('sample_labels', []) for i in self.headers]
        s = [[j.replace('_fwd', '') for j in i] for i in x]
        s = [[j.replace('_rev', '') for j in i] for i in s]
        return s


    def update_group_boundaries(self):
        """
        update group boundaries
        input: [0, 10, 30], [0, 5, 15]
        output: [0, 15, 45]
        """
        s = [i.get('group_boundaries', []) for i in self.headers]
        return reduce(lambda a,b: [i+j for i,j in zip(a,b)], s)


    def update_header(self):
        """
        update header line for final matrix file
        - sample_labels
        - group_boundaries
        """
        h = self.headers[0] # first header
        # update labels, boundaries
        s, g = self.check_labels()
        gb = self.update_group_boundaries()
        h.update({
            'sample_labels': s,
            'group_labels': g,
            'group_boundaries': gb
        })
        # return json.dumps(h)
        return h


    def run(self):
        if os.path.exists(self.o) and not self.overwrite:
            log.info('Matrix_rbind().run() skipped, file exists: {}'.format(self.o))
            return None
        # update header
        h = self.update_header()
        hs = json.dumps(h)
        # load matrix body, merge all
        mx = self.load_matrix2(self.m, header_only=False) # list of dict, key=group_labels
        with xopen(self.o, 'wt') as w:
            # save header
            w.write('@'+hs+'\n')
            # save matrix body
            for i in h.get('group_labels'):
                for j in mx:
                    w.write('\n'.join(j.get(i, []))+'\n')


def get_args():
    example = 'python matrix_rbind.py -m sens.gz anti.gz -o merge.gz'
    parser = argparse.ArgumentParser(
        prog='matrix_rbind.py',
        description='merge sens and anti matrix',
        epilog=example,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-m', dest='m', nargs='+', required=True,
        help='matrix files, by computeMatrix ')
    parser.add_argument('-o', dest='o', required=True,
        help='output matrix')
    parser.add_argument('-O', '--overwrite', dest='overwrite',
        action='store_true', help='Overwrite output file')
    return parser


def main():
    args = vars(get_args().parse_args())
    Matrix_rbind(**args).run()


if __name__ == '__main__':
    main()

# EOF