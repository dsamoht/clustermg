#!/usr/bin/env python
"""
usage:

Write a contig annotation table (contig_annotation.tsv) with results from Pfam and/or Kegg Ortholog.
"""
import pandas as pd
import os.path
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-p", "--pfamPath", required=False)
ap.add_argument("-k", "--keggPath", required=False)
ap.add_argument("-l", "--koList", required=False)
args = vars(ap.parse_args())

namesList = ["contigId", "accession_target", "tlen",
                    "query_name", "accession_query", "qlen",
                    "full_E_value", "full_score", "full_bias",
                    "dom_#", "dom_of", "dom_c_Evalue",
                    "dom_i_Evalue", "dom_score", "dom_bias",
                    "hmm_coord_from", "hmm_coord_to",
                    "ali_coord_from", "ali_coord_to",
                    "env_coord_from", "env_coord_to", "acc"]

# Pfam
if os.path.isfile(args['pfamPath']):
    df_pfam = pd.read_csv(args['pfamPath'], sep=" ", skiprows=3, skipfooter=10, header=None,
                    skipinitialspace=True, usecols=list(range(0, 22)),
                    names=namesList, na_values="-", engine="python")
    df_pfam = df_pfam.sort_values(by="full_E_value")
    df_pfam = df_pfam.drop_duplicates(subset=["contigId"], keep="first")
    df_pfam = df_pfam.iloc[:,[0,2,3,4,5,6]]
    df_pfam.columns = ["contigId", "tlen", "pfam_name", "pfam_id", "pfam_qlen", "pfam_E_value"]
else:
    df_pfam = pd.DataFrame(columns=["contigId", "tlen", "pfam_name", "pfam_id", "pfam_qlen", "pfam_E_value"])


# Kegg
if os.path.isfile(args['keggPath']):
    def swap_columns(df, col1, col2):
        col_list = list(df.columns)
        x, y = col_list.index(col1), col_list.index(col2)
        col_list[y], col_list[x] = col_list[x], col_list[y]
        df = df[col_list]
        return df
    
    ko_list = pd.read_csv(args['koList'], sep="\t")
    df_kegg = pd.read_csv(args['keggPath'], sep=" ", skiprows=3, skipfooter=10, header=None,
                    skipinitialspace=True, usecols=list(range(0, 22)),
                    na_values="-", engine="python")
    df_kegg = swap_columns(df_kegg, 3, 4)
    df_kegg.columns = namesList
    df_kegg = pd.merge(left=df_kegg, right=ko_list.iloc[:,[0, 11]], left_on="accession_query", right_on="knum")
    df_kegg.query_name = df_kegg.definition
    df_kegg = df_kegg.drop(["knum", "definition"], axis=1)
    df_kegg = df_kegg.sort_values(by="full_E_value")
    df_kegg = df_kegg.drop_duplicates(subset=["contigId"], keep="first")
    df_kegg = df_kegg.iloc[:,[0,2,3,4,5,6]]
    df_kegg.columns = ["contigId", "tlen", "kegg_name", "kegg_id", "kegg_qlen", "kegg_E_value"]
else:
    df_kegg = pd.DataFrame(columns=["contigId", "tlen", "kegg_name", "kegg_id", "kegg_qlen", "kegg_E_value"])

# Merge
df_annot = pd.merge(left=df_pfam, right=df_kegg, left_on="contigId", right_on="contigId", how="outer")
df_annot = df_annot.rename(columns={"tlen_x":"tlen"})
df_annot = df_annot.drop("tlen_y", axis=1)
df_annot.to_csv("contig_annotation.tsv", sep="\t")
