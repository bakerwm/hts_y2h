#!/usr/bin/env python3 

"""
Annotation alignment by gene_name/gene_id
output:
col-1: read_id 
col-2: gene_id
col-3: gene_name
col-4: position (====., ===.=, ==.==, =.===, .====)

Input:
1. bed6
2. gene_bed (tabix indexed, bgzipped)

Output:
1. txt file
"""

import os
import sys
from xopen import xopen
# from hiseq.utils.tabix_utils import Tabix
from tabix_utils import Tabix  # pip install pytabix


def find_pos(s, gene_list):
    """
    Find the position within gene_list

    Position:
    1. 3end: >>>>., .<<<<
    2. 5end: .>>>>, <<<<.
    3. genebody: >>.>>, <<.<<
    """
    pos = None
    if len(gene_list) > 0:
        for i in gene_list: # list_in_list, nested
            # skipped 
            if isinstance(pos, list):
                break
            if isinstance(i, list) and len(i) >= 6: # bed6
                if not s[5] == i[5]:
                    continue # require same strand
                # within gene
                if s[1] >= i[1] and s[2] <= i[2]:
                    gsize = int(i[2]) - int(i[1])
                    r1 = (int(s[1]) - int(i[1])) / gsize # left
                    r2 = (int(i[2]) - int(s[2])) / gsize # right
                    if r1 < 0.25:
                        p = "3end" if s[5] == "-" else "5end"
                    elif r2 < 0.25:
                        p = "5end" if s[5] == "-" else "3end"
                    else:
                        p = "genebody"
                    pos = [p, i[3]] # pos, gene_name
                # in IGR
                # pass
    return pos # 


def bed_to_gene(gene_bed, query_bed, anno_bed, overwrite=False):
    """
    Return the nearest gene to the query_bed

    Parameters:
      gene_bed (str): path to BED6 file, bgzipped, Tabix indexed; for genes
      query_bed (str): path to BED6 file, bgzipped, Tabix, for query
      anno_bed (str): output file, with pas updated
    """
    if os.path.exists(anno_bed) and overwrite is False:
        print('file exists, anno_bed() skipped, {}'.format(anno_bed))
        return None
    # iterate all genes
    # multiple-threads, if too-much bed (>1000000) 
    k = 0 # skipped lines
    gene_ix = Tabix(gene_bed) # gene idxed
    with xopen(query_bed) as r, xopen(anno_bed, 'wt') as w: # bed6 
        for l in r:
            # remove chr 
            l = l.replace('chr', '')
            s = l.strip().split() # query bed6 
            if len(s) < 6:
                k += 1 # record skipped lines
                continue # skipped
            # extract PAS within gene
            q = f'{s[0]}:{s[1]}-{s[2]}' # chr1:1-100
            # print('!A-1', q)
            gene_list = gene_ix.fetch(q)
            if gene_list is None:
                pos = ['null', 'null']
            else:
                gene_list = list(gene_list) # convert
                pos = find_pos(s, gene_list) # [gene_name, pos]
                if pos is None:
                    pos = ['null', 'null']
            w.write('\t'.join(s+pos)+'\n')


def main():
    if len(sys.argv) < 4:
        print("Usage: python anno_bed.py <gene.bed.gz> <query.bed> <out.bed>")
        sys.exit(1)
    gene_bed, query_bed, anno_bed = sys.argv[1:4]
    bed_to_gene(gene_bed, query_bed, anno_bed, True)


if __name__ == '__main__':
    main()

