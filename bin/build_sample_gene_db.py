#!/usr/bin/env python

import argparse
from pathlib import Path
import pandas as pd


class GeneInfoTable:

    def __init__(self, name: str, faa: str | Path, fcounts: str | Path):
        self.name: str = name
        self.faa: str | Path = Path(faa)
        self.fcounts: str | Path = Path(fcounts)
        self.info: list = []
        self.info_df: pd.DataFrame = None
        self.parse_faa()
        self.parse_fcounts()

    @staticmethod
    def fasta_generator(fasta_file: str | Path):
        with open(fasta_file, 'r') as fasta_handle:
            header = None
            seq = []
            for line in fasta_handle:
                if line.startswith('>'):
                    if header:
                        yield header, ''.join(seq)
                    header = line.strip().lstrip('>')
                    seq = []
                else:
                    seq.append(line.strip())
            if header:
                yield header, ''.join(seq)

    def parse_faa(self) -> None:
        self.info = [{
            'id': f"{self.name}_{i + 1}",
            'source': self.name,
            'header': header,
            'prodigal_id': header.split("ID=")[1].split(";")[0],
            'seq': seq,
            'contig': '_'.join(header.split()[0].split('_')[:-1]),
            'cds_len': 'NA',
            'mapped_reads': 'NA'
        } for i, (header, seq) in enumerate(self.fasta_generator(self.faa))]
        self.info_df = pd.DataFrame(self.info)

    def parse_fcounts(self) -> None:
        fcounts_df = pd.read_csv(self.fcounts, sep='\t', skiprows=1)
        fcounts_dict = fcounts_df.set_index('Geneid')[['Length', fcounts_df.columns[-1]]].to_dict(orient='index')
        self.info_df[['cds_len', 'mapped_reads']] = self.info_df['prodigal_id'].map(
            lambda x: fcounts_dict.get(x, ('NA', 'NA'))
        ).apply(pd.Series)

    def save_df(self) -> None:
        self.info_df.drop(columns=['prodigal_id'], inplace=True)
        self.info_df.to_csv(f'{self.name}_gene_info.tsv', sep='\t', index=False)


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', type=str, help='sample name', required=True)
    parser.add_argument('--faa', type=str, help='input faa file', required=True)
    parser.add_argument('--fcounts', type=str, help='input featureCounts file', required=True)
    args = parser.parse_args()
    return args


def main():
    args = arg_parser()
    name = args.name
    faa_file = args.faa
    fcounts_file = args.fcounts

    GeneInfoTable(name, faa_file, fcounts_file).save_df()


if __name__ == '__main__':
    main()