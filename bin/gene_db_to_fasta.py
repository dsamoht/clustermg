#!/usr/bin/env python

import argparse
from pathlib import Path

import pandas as pd


def make_faa(df: pd.DataFrame, name) -> None:
    with open(f'{name}.faa', 'w') as faa_out:
        for _, row in df.iterrows():
            faa_out.write(f'>{row['id']}\n{row["seq"]}\n')

def read_db(db_file: str | Path) -> pd.DataFrame:
    df = pd.read_csv(db_file, sep='\t', header=0)
    return df

def make_db_faa(db_df) -> None:
    df = db_df.loc[db_df['contig'].isna()]
    make_faa(df, 'db')

def make_sample_faa(db_df) -> None:
    df = db_df.loc[~db_df['contig'].isna()]
    make_faa(df, 'sample')

def arg_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--db', type=str, help='gene database', required=True)
    args = parser.parse_args()
    return args

def main():

    database = Path(arg_parser().db)
    df = read_db(database)
    make_db_faa(df)
    make_sample_faa(df)



if __name__ == "__main__":
    main()
