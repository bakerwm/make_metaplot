#!/usr/bin/env python
#-*- encoding: utf-8 -*-
"""
Download TSS from ENSEMBL BioMart, release-102
convert to BED format

How to?
1. remove chromosomes, `CHR_*` 
2. add "chr" to chromosome name
3. rename gene_name: "gene_name:chr:TSS"
4. convert to BED6 format, tss-start, tss-end

## see script: fetch_ensembl_tss.py
# Download TSS records from Ensembl-BioMart #
1. access to: http://nov2020.archive.ensembl.org/biomart/martview
2. set parameters: 
    CHOOSE DATABASE: "Ensembl Genes 102"
    CHOOSE DATASET: "Mouse genes (GRCm38.p6)"
3. choose "Attributes" on left-panel in the following order:
    "Features -> GENE"
    - Chromosome/scaffold name
    - Gene start (bp) 
    - Gene end (bp) 
    - Gene stable ID
    - Gene name
    - Strand
    - Transcription start site (TSS)
4. click "Count" (on top of left panel)
    update, Dataset 56305/56305 Genes 
5. Click "Results" next to "Count" button 
    set parameters on right panel,
    Export all results to "File", "CSV"
    click the option: "Unique results only"
6. Click "Go" button, the results will save to file "mart_export.txt"

> Here are the top-2 lines of the file:
Chromosome/scaffold name,Gene start (bp),Gene end (bp),Gene stable ID,Gene name,Strand,Transcription start site (TSS)
MT,15356,15422,ENSMUSG00000064372,mt-Tp,-1,15422
MT,15289,15355,ENSMUSG00000064371,mt-Tt,1,15289

> Here are the top-2 lines of output file
chrMT   15421   15422   ENSMUSG00000064372      mt-Tp:chrMT:15422       -
chrMT   15288   15289   ENSMUSG00000064371      mt-Tt:chrMT:15289       +

> Number of TSS records:
133,379 Mouse GRCm38.p6 Ensembl_release-102
206,880 Human GRCh38.p13 Ensembl_release-102
"""


import os
import sys
# from xopen import xopen


def format_tss(x, y):
    """
    How to?
    1. remove chromosomes, `CHR_*` 
    2. add "chr" to chromosome name
    3. rename gene_name: "gene_name:chr:TSS"
    4. convert TSS records to BED6 format, tss-start, tss-end
    
    Parameters:
    -----------
    x : str
        file, TSS records in csv format, download from Ensembl-BioMart
        see manual at bottom lines.
    y : str
        file, save the TSS BED file

    # Download TSS records from Ensembl-BioMart #
    1. access to: http://nov2020.archive.ensembl.org/biomart/martview
    2. set parameters: 
        CHOOSE DATABASE: "Ensembl Genes 102"
        CHOOSE DATASET: "Mouse genes (GRCm38.p6)"
    3. choose "Attributes" on left-panel in the following order:
        "Features -> GENE"
        - Chromosome/scaffold name
        - Gene start (bp) 
        - Gene end (bp) 
        - Gene stable ID
        - Gene name
        - Strand
        - Transcription start site (TSS)
    4. click "Count" (on top of left panel)
        update, Dataset 56305/56305 Genes 
    5. Click "Results" next to "Count" button 
        set parameters on right panel,
        Export all results to "File", "CSV"
        click the option: "Unique results only"
    6. Click "Go" button, the results will save to file "mart_export.txt"

    > Here are the top-2 lines of the file:
    Chromosome/scaffold name,Gene start (bp),Gene end (bp),Gene stable ID,Gene name,Strand,Transcription start site (TSS)
    MT,15356,15422,ENSMUSG00000064372,mt-Tp,-1,15422
    MT,15289,15355,ENSMUSG00000064371,mt-Tt,1,15289

    > Here are the top-2 lines of output file
    chrMT   15421   15422   ENSMUSG00000064372      mt-Tp:chrMT:15422       -
    chrMT   15288   15289   ENSMUSG00000064371      mt-Tt:chrMT:15289       +

    > Number of TSS records:
    133,379 Mouse GRCm38.p6 Ensembl_release-102
    206,880 Human GRCh38.p13 Ensembl_release-102
    """
    # duplicated names not allowed
    d = {}
    i = 0
    with open(x) as r, open(y, 'wt') as w:
        for l in r:
            i += 1
            if l.startswith('Chromosome') or l.startswith('CHR_'):
                continue
            p = l.strip().split(',') # csv file
            p[0] = 'chr'+p[0] # add chr to chromosome name
            tss = int(p[6]) # TSS in column-7
            p[1] = str(tss-1) # TSS start
            p[2] = str(tss) # TSS end
            p[3] = '{}:{}:{}'.format(p[4],p[0],p[6]) # name
            p[5] = '+' if p[5] == '1' else '-' if p[5] == '-1' else '*'
            b = p[:6] # BED6
            # check duplicated names
            dn = d.get(p[3], '')
            if dn == p[3]:
                print('line-{} found duplicated name: {}'.format(i, dn))
                continue # skip this line
            d.update({p[3]:p[3]}) # save id
            if y.endswith('.gtf'):
                b = bed2gtf(b, feature='gene')
            w.write('\t'.join(b)+'\n')


def main():
    if len(sys.argv) < 3:
        print('Usage: get_tss.py tss.csv tss.bed')
        sys.exit(1)
    format_tss(sys.argv[1], sys.argv[2])


if __name__ == '__main__':
    main()

# EOF