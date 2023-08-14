#!/usr/bin/env python 
#-*- encoding: UTF-8 -*-
"""
Parse genome annotation file (GTF), Ensembl, Gencode
optional: BED6, BED12 or BED15 format?

Filtering
1. by feature (gene)
2. by flanking distance
3. by biotype (protein_coding) 
4. by gene_size 

# for Ensembl-release 102; Gencode release-M25
# 110 genes were duplicated
# 
# keep the smallest gene_version:
Ptp4a1 Arhgef4 Gm24826 Gm26457 Septin2 Zfp813-ps Gm23370 Gm28040 Zc3h11a 
Gm28724 Snord80 Gm16701 Ndor1 Ndor1 Snora43 Nron Gm24350 Olfr1073-ps1 
Gm22813 A530058N18Rik Gm4430 Nkx2-2os Gm23925 1600017P15Rik Gm7270 Gm20690
Gm25820 4933434E20Rik Gm18433 Terc Tmigd3 Gm13301 Rmrp Gm25053 Pakap Pakap
Gm26379 Gm26047 Snora16a Gm22897 Gm28710 Jakmip1 Gm27680 Fam220a Gm24105 
Dlx6os1 Gm16499 Rprl1 Gm27013 Tmem147os Lim2 Mir1839 Olfr290 Aldoa Gm23128
Gm23604 Gm18645 Gm26265 Dpep2 Gm23377 Gm38642 Zkscan7 Gm16364 Rnu3a
C730027H18Rik Gm6729 Ddit3 Rnu3b1 Rnu3b3 Rnu3b4 Rnu3b2 Gm22711 Gm26413 
St6galnac2 1700030C10Rik Gm22149 Gm25203 Gm27528 Gm35558 Gm35558 Gm24022 
Ighv5-8 Ighv1-13 Vmn1r216 Gm36638 Gm9025 Nnt Vmn2r-ps111 Gm6740 Gm5089 
4930594M22Rik Mirt2 Gcat Gm27825 Gm28023 Gm41392 Gm38619 Gm25617 Atp5o Gm29719
Gm5966 Snhg4 Mir1949 Pcdha11 Arhgap26 Gm23927 Zfp91 Gm35438 Gm17522 Gm23786

# number of genes
source           ensembl    gencode
release          102        M25
date             2020-11    2020-03-24
genes            55401      55401
protein_coding   21859      21859

filtered:
genes            55367      55291 (remove duplicated gene_names)
protein_coding   21898      21831 (only protein_coding)
"""


import os
import sys
import io
import pathlib
import argparse
import pybedtools
from datetime import datetime
from xopen import xopen


def cur_time():
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def filt_gtf(gtf, bed, **kwargs):
    """
    Filt Ensembl GTF file
    - feature: gene|exon|CDS|transcript, default: None
    - gene_biotype: TEC|protein_coding, default: None 
    - gene_size: --
    - flanking_size: --
    duplicated genes: Zfp91, Zkscan7, 
    """
    overwrite = kwargs.get('overwrite', False)
    append_biotype = kwargs.get('append_biotype', False)
    if os.path.exists(bed) and overwrite is False:
        print('file exists, gtf_to_bed() skipped, {}'.format(bed))
        return None
    # check duplicated genes
    g = {}
    for i in read_gtf(gtf, **kwargs):
        # bed6 + gene_biotype + gene_version
        pre = g.get(i[3], None)
        if pre is None:
            g.update({i[3]:i})
        else:
            # compare gene_version
            # print('!A-1', i[3])
            if i[-1] > pre[-1]: #!!! keep the largest version !!!
                g.update({i[3]:i})
    # save
    with xopen(bed, 'wt') as w:
        for k,v in g.items():
            # update bed-start
            v[1] = str(int(v[1]) - 1)
            bed6 = v[:-1] if append_biotype else v[:-2]            
            w.write('\t'.join(bed6)+'\n')


def read_gtf(gtf, **kwargs):
    """
    Read Ensembl GTF file

    Parameters:
    -----------
    feature : str
        gene|exon|CDS|transcript, default: None
    gene_biotype : str
        TEC|protein_coding, default: None
    Convert GTF to BED format
    filt by: gene_biotype (ENSEMBL), gene_type (GENCODE)

    duplicate gene records found in Ensembl.
    eg: Zfp91, Zkscan7, ...
    filtered by: "gene_version" !!!
    """
    overwrite = kwargs.get('overwrite', False)
    line = 0
    with xopen(gtf) as r:
        for g in parse_gtf(r, **kwargs):
            line += 1
            # add chr to chromosome name
            chr = g[0] # if g[0].startswith('chr') else 'chr'+g[0]
            n = parse_gtf_desc(g[8], 'gene_name') # column-9
            # gene version in "GenCODE"
            i = parse_gtf_desc(g[8], 'gene_id')
            if '.' in i: # gencode
                i,y = i.split('.', 1)
                v = int(y)
            else: # Ensembl
                v = parse_gtf_desc(g[8], 'gene_version')
                v = 1 if v is None else int(v)
            n = i if n is None else n
            # warning: if n is None
            if n is None:
                print('line-{} failed, could not found'.format(line))
            bt = parse_gtf_desc(g[8], 'gene_biotype')
            # 6-columns + gene_biotype + gene_version
            yield [chr, g[3], g[4], n, '255', g[6], bt, v]


def parse_gtf(fh, **kwargs):
    """
    filt GTF records by [feature] and [gene_biotype]

    Parameters:
    -----------
    feature : str
        gene|exon|CDS|transcript, default: None (all)
    gene_biotype : str
        TEC|protein_coding, default: None (all)
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
    read Description column of GTF (column-9), gene_id "ENSMUSG00000118491"; ...
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


def filt_by_size(x, y, size=5000, overwrite=False):
    """
    gene_biotype saved in column-7 (last column)
    """
    if os.path.exists(y) and overwrite is False:
        print('file exists, flit_bed_size() skipped, {}'.format(y))
    else:
        bed = pybedtools.BedTool(x).filter(lambda i: len(i) >= size) # gene size
        bed.saveas(y)


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
            if flank_size < 0:
                pass # save all
            elif i.strand == '+':
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


def filt_by_biotype(x, y, gene_biotype=None, overwrite=False):
    """
    gene_biotype saved in column-7 (last column)
    """
    if os.path.exists(y) and overwrite is False:
        print('file exists, flit_bed_size() skipped, {}'.format(y))
    else:
        bed = pybedtools.BedTool(x)
        if isinstance(gene_biotype, str):
            bed = bed.filter(lambda i: i[-1] == gene_biotype)
        elif isinstance(gene_biotype, list):
            bed = bed.filter(lambda i: i[-1] in gene_biotype)
        else:
            pass
        bed.saveas(y)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='gtf', required=True, 
        help='The gene GTF file')
    parser.add_argument('-o', dest='out_dir', required=False, 
        help='The out_dir')
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
    2. genebody size (750)
    3. flanking size 
    4. protein coding
    """
    args = get_args().parse_args()
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
    bed_pc = os.path.join(out_dir, prefix+'.gene.f{}.pc.bed'.format(fs))
    bed_size = os.path.join(out_dir, prefix+'.gene.f{}.pc.g{}.bed'.format(fs, gs))

    # step1. gtf to bed, feature, (all genes)
    print('[{}] 1. extract all gene'.format(cur_time()))
    kwargs = {
        'feature': args.feature, 
        'gene_biotype': None, # args.gene_biotype,
        'append_biotype': True, # args.append_biotype,
        'overwrite': args.overwrite,
        }
    filt_gtf(args.gtf, bed_gene, **kwargs) # protein_coding, genes

    # step2. filt by flanking size
    print('[{}] 2. filter by, flank_size={}'.format(cur_time(), flank_size))
    filt_by_flank(bed_gene, bed_flank, flank_size=flank_size, overwrite=args.overwrite)

    # step3. filt by genetype (protein_coding)
    print('[{}] 3. filter by, gene_biotype={}'.format(cur_time(), args.gene_biotype))
    filt_by_biotype(bed_flank, bed_pc, gene_biotype=args.gene_biotype, overwrite=args.overwrite)

    # step4. filt by gene size
    print('[{}] 4. filter by, gene_size={}'.format(cur_time(), gene_size))
    filt_by_size(bed_pc, bed_size, size=args.gene_size, overwrite=args.overwrite)

    print('Finish, save to file: {}'.format(bed_size))

if __name__ == '__main__':
    main()
    
# EOF