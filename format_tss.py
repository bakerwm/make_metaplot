#!/usr/bin/env python
"""
Download TSS from ENSEMBL BioMart, release-102
save as BED
column-1: chromosome name, add "chr"
column-2: gene start (1-index to 0-index)
column-3: gene end 
column-4: gene_id (ENSEMBL_id)
column-5: gene_name (SYMBOL) 
column-6: strand (1,-1 to +,-) 
column-7: TSS 
"""

import os
import sys

def main():
    if len(sys.argv) < 2:
        print('Usage: get_tss.py tss.csv > tss.bed')
        sys.exit(1)
    with open(sys.argv[1]) as r:
        for l in r:
            if l.startswith('Chromosome'):
                continue
            p = l.strip().split(',')
            p[0] = 'chr'+p[0] # update chr
            p[1] = str(int(p[1]) - 1) # to 0-index, BED-start
            p[5] = '+' if p[5] == '1' else '-' # fix strand
            print('\t'.join(p))


if __name__ == '__main__':
    main()

# EOF