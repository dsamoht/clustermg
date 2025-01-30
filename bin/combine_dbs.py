#!/usr/bin/env python

import argparse
import csv
from pathlib import Path

import pandas as pd


def combine_df(dbs: str | Path) -> None:
    db_to_concat = []
    for db in dbs:
        df = pd.read_csv(db, sep='\t', header=0)
        db_to_concat.append(df)
    combined_df = pd.concat(db_to_concat, ignore_index=True)
    combined_df.to_csv('gene_db.tsv', sep='\t', index=False, float_format='%.0f', quoting=csv.QUOTE_NONE)

def arg_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--dbs', type=str, help='gene databases', required=True, nargs='*')
    args = parser.parse_args()
    return args

def main():

    databases = arg_parser().dbs
    combine_df(databases)

if __name__ == "__main__":
    main()
