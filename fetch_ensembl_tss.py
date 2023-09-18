#!/usr/bin/env python3
#-*- encoding: utf8 -*-

"""
Name: fetch TSS from ensembl

Parameters:

release: 
organism:

in this example:
release-102: http://nov2020.archive.ensembl.org
database: ENSEMBL_MART_ENSEMBL
dataset:  mmusculus_gene_ensembl

How to?
1. remove chromosomes, `CHR_*` 
2. rename gene_name: "gene_name:chr:TSS"
3. convert TSS records to BED6 format, tss-start, tss-end
4. use `chromToUcsc` to convert files into UCSC annotation

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
    - Gene name
    - Gene stable ID
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
import biomart
from biomart import BiomartServer


def list_ensembl_arvhives():
    """
    return:
    name, date, url, version
    """
    pass


def list_ensembl():
    pass


def list_ensembl_genomes():
    """
    return:
    biomart, version
    """
    pass


def build_biomart(database, dataset, url):
    """
    helper functions:
    is_valid_database:
    is_valid_dataset:
    is_valid_organism:
    
    url: http://nov2020.archive.ensembl.org # release-102
    database: ENSEMBL_MART_ENSEMBL
    dataset:  mmusculus_gene_ensembl 
    """
    try:
        server = BiomartServer(url+'/biomart')
    except:
        print('Could not access Ensembl server')
        server = None
    if isinstance(server, biomart.BiomartServer):
        if server.is_alive:
            ensembl = server.databases[database]
            mart = ensembl.datasets[dataset]
            return mart


def download_tss(out_txt, oganism=None, release=None):
    """
    release-102: http://nov2020.archive.ensembl.org
    database: ENSEMBL_MART_ENSEMBL
    dataset:  mmusculus_gene_ensembl 
    """
    if os.path.exists(out_txt):
        print(f'file exists: {out_txt}')
        return out_txt
    mart = build_biomart(
        database = 'ENSEMBL_MART_ENSEMBL',
        dataset = 'mmusculus_gene_ensembl',
        url = 'http://nov2020.archive.ensembl.org'
    )
    # query, for TSS
    query = {
        'attributes': [
            'chromosome_name', 
            'start_position', 
            'end_position', 
            'external_gene_name',
            'ensembl_gene_id', 
            'strand',
            'transcription_start_site',
        ]
    }
    # search
    tss_n = mart.count(query)
    if tss_n > 0:
        tss = mart.search(query)
        out_tmp = out_txt + '.tmp'
        with open(out_tmp, 'wt') as w:
            w.write(tss.text)
        fix_tss(out_tmp, out_txt)
        if os.path.exists(out_tmp):
            os.remove(out_tmp)
    print(f'found {tss_n} TSS for organism')


def fix_tss(in_bed, out_bed):
    """
    1. remove chromosomes, `CHR_*` 
    2. rename gene_name: "gene_name:chr:TSS"
    3. convert TSS records to BED6 format, tss-start, tss-end
    """
    with open(in_bed) as r, open(out_bed, 'wt') as w:
        for l in r:
            if l.startswith('CHR_'):
                continue
            p = l.strip().split()
            p[3] = f'{p[3]}:{p[0]}:{p[6]}' # gene_name:chr:tss
            p[1] = str(int(p[6]) - 1)
            p[2] = p[6]
            p[5] = '-' if p[5] == -1 else '+'
            w.write('\t'.join(p[:6])+'\n')


def main():
    if len(sys.argv) < 3:
        print('Usage: python fetch_ensembl_tss.py <genome> <tss.txt>')
        sys.exit(1)
    genome, out_txt = sys.argv[1:3]
    download_tss(out_txt)


if __name__ == '__main__':
    main()
