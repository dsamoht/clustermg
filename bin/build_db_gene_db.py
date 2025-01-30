#!/usr/bin/env python

import argparse
from collections import defaultdict
from pathlib import Path
import pandas as pd


class GeneInfoTable:

    def __init__(self, name: str, faa: str | Path):
        self.name: str = name
        self.faa: str | Path = Path(faa)
        self.info : list = []
        self.info_df: pd.DataFrame = None
        self.parse_faa()

    @staticmethod
    def fasta_to_dict(fasta_file: str | Path):
        with open(fasta_file, 'r') as fasta_handle:
            header = None
            for line in fasta_handle:
                if line.startswith('>'):
                    if header:
                        yield header, ''.join(seq)
                    header = line.strip().split('>')[1]
                    seq = []
                else:
                    seq.append(line.strip())
            if header:
                yield header, ''.join(seq)

    def parse_faa(self) -> None:
        records = dict(self.fasta_to_dict(self.faa))
        self.info = [{
            'id': f"{self.name}_{i+1}",
            'source': self.name,
            'header': header,
            'seq': seq
        } for i, (header, seq) in enumerate(records.items())]
        self.info_df = pd.DataFrame(self.info)
    
    def save_df(self) -> None:
        self.info_df.to_csv(f'{self.name}_gene_info.tsv', sep='\t', index=False)

def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', type=str, help='database name', required=True)
    parser.add_argument('--faa', type=str, help='input faa file', required=True)
    args = parser.parse_args()
    return args

def main():
    args = arg_parser()
    name = args.name
    faa_file = args.faa

    GeneInfoTable(name, faa_file).save_df()

if __name__ == '__main__':
    main()