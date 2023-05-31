#!/usr/bin/env python 
#-*- encoding: UTF-8 -*-
"""
Update genes by PAS 

Optional:
1. distal PAS within genebody (in 3' UTR) 
2. distal PAS within termination zone (x kb downstream of PAS/TES)
3. top PAS-seq signal (not availabel)
"""


import os
import sys
from xopen import xopen
from hiseq.utils.tabix_utils import Tabix


def distal_pas_in_tz(x, pas, y, overwrite=False):
    """
    Distal PAS in termination zone, 3' flanking region
    """
    pass


def pick_distal_pas(gene, pas_list):
    """
    Parameters:
      gene (list): bed6 
      pas_list (list): list of list for PAS annotation
    """
    pas = None
    if len(pas_list) > 0:
        for i in pas_list: # list_in_list (nested)
            if isinstance(i, list) and len(i) >= 6:
                if pas is None:
                    pas = i[2] if i[5] == '+' else i[1] # first record
                if i[5] == '-' and i:
                    pas = i[1] if pas > i[1] else pas
                else:
                    pas = i[2] if pas < i[2] else pas
    return pas
                


def distal_pas_in_gb(gene_bed, pas_bed, gene_pas, tes_for_na=True, overwrite=False):
    """
    Distal PAS in genebody, 3' UTR

    Parameters:
      gene_bed (str): path to BED6 file, bgzipped, Tabix; for genes
      pas_bed (str): path to BED6 file, bgzipped, Tabix, for PAS annotation 
      gene_pas (str): output file, with pas updated 
      tes_for_na (bool): use TES instead, if no pas found for the gene (True), 
        or discard the no-pas genes
    """
    if os.path.exists(gene_pas) and overwrite is False:
        print('file exists, flit_by_pas() skipped, {}'.format(gene_pas))
        return None
    # iterate all genes
    # multiple-threads, if too-much bed (>1000000) 
    k = 0 # skipped lines
    pas_ix = Tabix(pas_bed) #
    with xopen(gene_bed) as r, xopen(gene_pas, 'wt') as w: # bed6 
        for l in r:
            s = l.strip().split() # bed6 
            if len(s) < 6:
                k += 1
                continue # skipped
            # extract PAS within gene
            region = f'{s[0]}:{s[1]}-{s[2]}'
            try:
                pas_list = list(pas_ix.fetch(region)) # candidates
                pas = pick_distal_pas(s, pas_list)
            except:
                pas = None
                # print(l.strip())
            if not tes_for_na and pas is None:
                continue # skipped
            # update PAS
            if s[5] == '-':
                s[1] = s[1] if pas is None else pas
            else:
                s[2] = s[2] if pas is None else pas
            w.write('\t'.join(s)+'\n')
            

def main():
    if len(sys.argv) < 4:
        print(
        """
        Usage: python update_by_pas.py <gene.bed.gz> <pas.bed.gz> <out.bed>

        Parameters:
          gene.bed.gz   str    path to gene file, bed6 format
          pas.bed.gz    str    path to PAS file, bgzipped, tabix indexed, bed6 
          out.bed       str    path to output file, in BED format as gene.bed
        """
        )
        sys.exit(1)
    gene_bed, pas_bed, gene_pas = sys.argv[1:4]
    distal_pas_in_gb(gene_bed, pas_bed, gene_pas, tes_for_na=False, overwrite=True)


if __name__ == '__main__':
    main()
