#!/usr/bin/env python

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


class GeneClusterMatrix:

    def __init__(self, clstr: str | Path, gene_db: str | Path):
        self.clstr: str | Path = Path(clstr)
        self.gene_db: str | Path = Path(gene_db)
    
    def make_matrix(self):
        cluster_df = pd.read_csv(self.clstr, sep='\t', header=None, names=['representative', 'member'])
        gene_df = pd.read_csv(self.gene_db, sep='\t')
        total_mapped_reads_per_source = gene_df.groupby('source')['mapped_reads'].sum().to_dict()
        gene_df['RPKM'] = gene_df.apply(
            lambda row: (row['mapped_reads'] * 1e9) / (row['cds_len'] * total_mapped_reads_per_source[row['source']]),
            axis=1)
        member_to_rep = dict(zip(cluster_df['member'], cluster_df['representative']))
        rpkm_matrix = defaultdict(lambda: defaultdict(float))

        for _, row in gene_df.iterrows():
            gene_id = row['id']
            source = row['source']
            rpkm = row['RPKM']

            if gene_id in member_to_rep:
                representative = member_to_rep[gene_id]
                rpkm_matrix[representative][source] += rpkm

        rpkm_df = pd.DataFrame.from_dict(rpkm_matrix, orient='index')
        rpkm_df = rpkm_df.fillna(0)
        rpkm_df.to_csv('mmseqs_cluster_rpkm.tsv', sep='\t')


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--clstr', type=str, help='mmseqs tsv output', required=True)
    parser.add_argument('--genedb', type=str, help='gene database', required=True)
    args = parser.parse_args()
    return args


def main():
    args = arg_parser()
    clstr = args.clstr
    genedb = args.genedb
    GeneClusterMatrix(clstr, genedb).make_matrix()


if __name__ == '__main__':
    main()
