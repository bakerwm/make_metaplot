#!/usr/bin/env python 
"""
Purpose: extract genes with APA site

Reference: 
1.Kamieniarz-Gdula, K. et al. Selective Roles of Vertebrate PCF11 in Premature and Full-Length Transcript Termination. Molecular Cell 74, 158-172.e9 (2019).

Criteria:
1. extract all gene
2. filt by PAS DB (distal PAS site), 3' mRNA-seq; polyA_DB, polyASite
3. distance to flanking gene (6 kb)
4. genebody size (5 kb)
5. filt by features (protein coding)

Hg19/GRCh37 was used as the reference genome. GENCODE release 19 was used for 
gene annotations: https://www.gencodegenes.org/. This annotation includes
57820 genes (20345 protein-coding, 37475 non-coding). For downstream analysis,
we selected a subset of 11947 genes (9095 protein-coding, 2852 non-coding) 
that satisfied all of the following 3 criteria: 1) had at least one active
PAS (see 30 mRNA-seq analysis below for details); 2) did not overlap with 
another annotated gene on the same strand; 3) had a 30 end isolated by at 
least 6 kb from the downstream annotated gene on the same strand. Those 
strand-specific isolation criteria allowed to unambiguously assign the 
directional RNA-seq signal (chrRNA-seq, mNET-seq and 30 mRNA-seq) to the 
end of each gene, and also to compute distal alternative polyadenylation 
(APA) downstream of annotated gene ends (see below). 6 kb isolation was
used because visual inspection of the data in genome browser revealed 
usage of cryptic non-annotated PASs used upon PCF11 depletion within 
this window. For meta-profiles and heatmaps, a subset of protein-coding
genes longer that 5 kb was used (n = 8389), or a further subset of those 
as indicated in the figure legend. For calculation of distances between 
genes (Figures 4 and 7) the downstream distance from the gene’s 30 
end to any other annotated gene end (50 or 30) on either strand was computed.
"""

import os
import sys
import io
import pathlib
import argparse
import pybedtools
from datetime import datetime
from xopen import xopen


def gtf_to_bed(gtf, bed, **kwargs):
    """
    Parameters:
    -----------
    feature : str
        gene|exon|CDS|transcript, default: None
    gene_biotype : str
        TEC|protein_coding, default: None
    Convert GTF to BED format
    filt by: gene_biotype (ENSEMBL), gene_type (GENCODE)
    """
    overwrite = kwargs.get('overwrite', False)
    if os.path.exists(bed) and overwrite is False:
        print('file exists, gtf_to_bed() skipped, {}'.format(bed))
    else:
        with xopen(gtf) as r, xopen(bed, 'wt') as w:
            for g in read_gtf(r, **kwargs):
                g_chr = g[0] if g[0].startswith('chr') else 'chr'+g[0]
                g_start = str(int(g[3])-1)
                g_name = parse_gtf_desc(g[8], 'gene_name')
                if g_name is None:
                    g_name = parse_gtf_desc(g[8], 'gene_id') # alternative
                gene_biotype = parse_gtf_desc(g[8], 'gene_biotype')
                bed6 = [g_chr, g_start, g[4], g_name, '254', g[6]]
                if kwargs.get('append_biotype', False):
                    bed6.append(gene_biotype)
                bed6 = list(map(str, bed6))
                w.write('\t'.join(bed6)+'\n')


def read_gtf(fh, **kwargs):
    """
    Parameters:
    -----------
    feature : str
        gene|exon|CDS|transcript, default: None (all)
    gene_biotype : str
        TEC|protein_coding, default: None (all)

    filt GTF records by [feature] and [gene_biotype]
    """
    feature = kwargs.get('feature', None)
    gene_biotype = kwargs.get('gene_biotype', None)
    for l in fh:
        p = l.strip().split('\t')
        # filtering
        if len(p) < 9:
            continue
        # filt by feature
        if not feature is None:
            if isinstance(feature, str) and not p[2] == feature:
                continue
            elif isinstance(feature, list) and not p[2] in feature:
                continue
            else:
                pass
        # filt by gene_biotype
        gb = parse_gtf_desc(p[8], 'gene_biotype')
        if not gene_biotype is None:
            if isinstance(gene_biotype, str) and not gb == gene_biotype:
                continue
            elif isinstance(gene_biotype, list) and not gb in gene_biotype:
                continue
            else:
                pass
        yield p


def parse_gtf_desc(s, key='gene_name'):
    """
    Description in GTF:
    gene:
    ENSEMBL:
    gene_id "ENSMUSG00000102693"; gene_version "1"; gene_name "4933401J01Rik"; 
    gene_source "havana"; gene_biotype "TEC"; havana_gene "OTTMUSG00000049935";
    havana_gene_version "1";
    GENCODE:
    gene_id "ENSG00000227232.4"; transcript_id "ENSG00000227232.4"; 
    gene_type "pseudogene"; gene_status "KNOWN"; gene_name "WASH7P"; 
    transcript_type "pseudogene"; transcript_status "KNOWN"; 
    transcript_name "WASH7P"; level 2; havana_gene "OTTHUMG00000000958.1";
    """
    d = {}
    if isinstance(s, str):
        for p in s.strip().split(';'):
            if " " in p:
                a, b = p.split()[:2]
                b = b.replace('"', '')
                d.update({a:b})
    # fix gene_type (GENCODE), gene_biotype (ENSEMBL)
    if key == 'gene_type' and 'gene_biotype' in s:
        key = 'gene_biotype'
    elif key == 'gene_biotype' and 'gene_type' in s:
        key = 'gene_type'
    return d.get(key, None)


def filt_by_flank(x, y, flank_size=6000, overwrite=False):
    """
    Filt the BED file by flanking size
    1. flanking BED (6kb) not including genes
    """
    if os.path.exists(y) and overwrite is False:
        print('file exists, flit_by_flank() skipped, {}'.format(y))
        return None
    # sort bed
    xb = pybedtools.BedTool(x).sort()
    # filter by flank size
    with xopen(y, 'wt') as w:
        pre_gene = [None, None] #fwd, rev
        for i in xb:
            if i.strand == '+':
                if pre_gene[0] is None:
                    pre_gene[0] = i
                else:
                    if pre_gene[0].chrom == i.chrom and i.start - pre_gene[0].end <= flank_size:
                        continue
                    else:
                        pre_gene[0] = i
            elif i.strand == '-':
                if pre_gene[1] is None:
                    pre_gene[1] = i
                else:
                    if pre_gene[1].chrom == i.chrom and i.start - pre_gene[1].end <= flank_size:
                        continue
                    else:
                        pre_gene[1] = i
            bed6 = list(map(str, i[:7])) # bed6+1
            w.write('\t'.join(bed6)+'\n')


def filt_by_pas(x, pas, y, overwrite=False):
    """
    Filter by another BED (pas)
    extract the distal pas site: tss-pas
    """
    if os.path.exists(y) and overwrite is False:
        print('file exists, flit_by_pas() skipped, {}'.format(y))
        return None
    # intersect a and b, by strand
    xp = pybedtools.BedTool(x).intersect(pas, s=True, wo=True)
    xp = xp.sort() #
    # choose distal 3' PAS site
    with xopen(y, 'wt') as w:
        pre_pas = None
        pre_bed = None
        for i in xp:
            i_pas = int(i[9]) # PAS site, bed7+bed6
            if pre_bed is None: # init
                pre_pas = i_pas
                pre_bed = i
            else:
                if i[3] == pre_bed[3]: # same gene, update pas_pos
                    if i[5] == '-' and pre_pas > i_pas: #
                        pre_pas = i_pas
                    elif i[5] == '+' and pre_pas < i_pas: #
                        pre_pas = i_pas
                    else:
                        pass
                else: # save record
                    bed6 = pre_bed[:7]
                    if bed6[5] == '-':
                        bed6[1] = pre_pas
                    else:
                        bed6[2] = pre_pas                        
                    bed6 = list(map(str, bed6))
                    # fix pre-bed                    
                    if int(bed6[1]) < int(bed6[2]):
                        w.write('\t'.join(bed6)+'\n')
                    pre_pas = i_pas
                    pre_bed = i                    
        # last record
        bed6 = pre_bed[:7]
        if bed6[5] == '-':
            bed6[1] = pre_pas
        else:
            bed6[2] = pre_pas
        bed6 = list(map(str, bed6))
        w.write('\t'.join(bed6)+'\n')


def filt_bed_size(x, y, gene_size=5000, gene_biotype=None, overwrite=False):
    """
    gene_biotype saved in column-7 (last column)
    """
    if os.path.exists(y) and overwrite is False:
        print('file exists, flit_bed_size() skipped, {}'.format(y))
    else:
        bed = pybedtools.BedTool(x).filter(lambda i: len(i) > gene_size) # gene size
        if isinstance(gene_biotype, str):
            bed = bed.filter(lambda i: i[-1] == gene_biotype)
        elif isinstance(gene_biotype, list):
            bed = bed.filter(lambda i: i[-1] in gene_biotype)
        else:
            pass
        bed.saveas(y)


def cur_time():
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='gtf', required=True, 
        help='The gene GTF file')
    parser.add_argument('-o', dest='out_dir', required=False, 
        help='The out_dir')
    parser.add_argument('-pas', dest='pas_bed', required=False,
        help='The PAS file in BED format')
    parser.add_argument('-f', '--feature', 
        help='[gene|exon|CDS|transcript], default [None] all')
    parser.add_argument('-g', '--gene-biotype', dest='gene_biotype',
        help='The gene_biotype, [protein_coding|...], default: [None] all')
    parser.add_argument('-x', '--append-biotype', dest='append_biotype', 
        action='store_true', help='Add gene_biotype to column-7 in BED')
    parser.add_argument('-gs', '--gene-size', dest='gene_size', type=int, 
        default=500, help='The minimum size of BED feature, default [0]')
    parser.add_argument('-fs', '--flank-size', dest='flank_size', type=int, 
        default=1, help='The minimum distance to neighbour gene, default [0]')
    parser.add_argument('-O', '--overwrite', action='store_true',
        help='Overwrite exists file')
    return parser


def main():
    """
    1. extract all gene
    2. distance to flanking gene (6kb)
    3. filt by PAS DB (distal PAS site)
    4. protein_coding, and genebody size (5kb)
    """
    args = get_args().parse_args()
    # gtf = args.gtf
    gene_size = args.gene_size
    flank_size = args.flank_size
    out_dir = args.out_dir
    # out_dir
    if not isinstance(out_dir, str):
        out_dir = pathlib.Path.cwd()
    out_dir = os.path.abspath(out_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # prepare files
    gs = '{:d}k'.format(int(gene_size/1000)) if gene_size > 1000 else gene_size
    fs = '{:d}k'.format(int(flank_size/1000)) if flank_size > 1000 else flank_size
    prefix = os.path.splitext(os.path.basename(args.gtf))[0]
    prefix = os.path.splitext(prefix)[0] if prefix.endswith('.gtf') else prefix
    bed_gene = os.path.join(out_dir, prefix+'.gene.bed')
    bed_flank = os.path.join(out_dir, prefix+'.gene.f{}.bed'.format(fs))
    bed_pas = os.path.join(out_dir, prefix+'.gene.f{}.pas.bed'.format(fs))
    bed_size = os.path.join(out_dir, prefix+'.gene.f{}.pas.g{}.pc.bed'.format(fs, gs))
    # bed_pc = os.path.join(out_dir, prefix+'.gene.pas.f{}.g{}.pc.bed'.format(fs, gs))

    # step1. gtf to bed, feature, (all genes)
    print('[{}] 1. extract all gene'.format(cur_time()))
    kwargs = {
        'feature': args.feature, 
        'gene_biotype': None, # args.gene_biotype,
        'append_biotype': True, # args.append_biotype,
        'overwrite': args.overwrite,
        }
    gtf_to_bed(args.gtf, bed_gene, **kwargs) # protein_coding, genes

    # step2. filt by flanking size
    print('[{}] 2. filter by flank size={}'.format(cur_time(), flank_size))
    filt_by_flank(bed_gene, bed_flank, flank_size=flank_size, overwrite=args.overwrite)

    # step3. filt by PAS DB
    print('[{}] 3. filter by PAS_DB'.format(cur_time()))
    if os.path.exists(args.pas_bed):
        filt_by_pas(bed_flank, args.pas_bed, bed_pas, overwrite=args.overwrite)

    # step4. filt by gene size, gene_biotype=protein-coding
    print('[{}] 4. filter by, protein-coding and gene size={}'.format(cur_time(), gene_size))
    filt_bed_size(bed_pas, bed_size, gene_biotype='protein_coding', overwrite=args.overwrite)

    # print('Finish, save to file: {}'.format(bed_pas))

if __name__ == '__main__':
    main()
    
# EOF