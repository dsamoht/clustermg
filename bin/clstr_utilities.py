#!/usr/bin/env python
"""Utilities."""
from __future__ import annotations
from argparse import ArgumentParser
from pathlib import Path
import pickle
import polars as pl


class ClusterSeq:
    """
    Main class of clustering utilities.
    """
    def __init__(self, clstr: str | Path, fasta: str | Path, output_dir: str | Path) -> None:
        """
        An object containing:
        - DataFrame `seq_info_df` where columns are:
        `seq_id`: str, `sequence`: str, `cluster`: str, `source`: str, `representative`: boolean;
        - utility functions to filter `seq_info_df` and output files of different format.
        Parameters:
        -----------
        `clstr`: .clstr file (`CD-HIT` output);
        `fasta`: fasta file given as input to `CD-HIT`.
        """

        self.clstr_file = Path(clstr)
        self.clstr_base_name = self.clstr_file.name.rstrip(r".clstr")
        self.fasta_file = Path(fasta)
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        self.output_dir = Path(output_dir)
        self._sequence_clstr_info_df()
        self._matrix_from_cdhit()

    def _sequence_clstr_info_df(self) -> None:
        """
        Build `seq_info_df`.
        """
        sequence_info = {}
        with open(self.fasta_file, "r", encoding="UTF-8") as fasta_input:
            for line in fasta_input:
                if line.startswith(">"):
                    current_name = line.strip().split()[0]
                    sequence_info[current_name] = [line.strip(), ""]
                else:
                    sequence_info[current_name][1] += line.strip()

        clstr_info = {}
        with open(self.clstr_file, "r", encoding="UTF-8") as clstr_input:
            index = -1
            for line in clstr_input:
                rep = 'False'
                if line.startswith(">"):
                    current_clstr = f"cluster_{line.split()[1]}"
                else:
                    index += 1
                    if line.strip().endswith(r"*"):
                        rep = 'True'
                    seq_name, seq = sequence_info[line.split()[2].rstrip(r"...")]
                    clstr_info[index] = seq_name, seq, seq_name.split("|")[0][1:], current_clstr, rep

        clstr_info = {str(k):v for k,v in clstr_info.items()}
        self.seq_info_df = pl.from_dict(clstr_info)
        self.seq_info_df = self.seq_info_df.transpose(column_names=["seq_id", "seq", "source", "cluster", "representative"])
        self.seq_info_df = self.seq_info_df.with_columns(pl.col('representative')=='True')
        self.seq_info_df = self.seq_info_df.rename({"seq_id": "seq_name_fasta"})
        self.seq_info_df = self.seq_info_df.insert_at_idx(0, self.seq_info_df['seq_name_fasta'].str.extract(r"^>(\S*)", group_index=1).alias('seq_id'))
        contig_names = self.seq_info_df['seq_id'].str.replace(pattern=r"(_[0-9]*$)", value='')
        contig_names = contig_names.rename("contigId")
        self.seq_info_df = self.seq_info_df.insert_at_idx(1, contig_names)
        output_file = self.output_dir.joinpath(f"clusters_info.tsv")
        self.seq_info_df.write_csv(output_file, separator='\t')


    def contiguity_network(self, file: str | Path) -> None:
        """
        Write a contiguity network to `file`, where nodes are clusters number
        and edges are contiguous genes.

        How to interpet:
        If clusters A and B are joined by 1 edge, it means that 1 gene from
        cluster A is contiguous to 1 gene from cluster B. A chain of clusters
        might represent a genomic region of interest
        (example: biosynthetic gene clusters).

        Parameter:
        -----------
        `file`: path to output file
        """
        tmp_df = self.seq_info_df.filter(pl.col('source')=='SAMPLE')



    def clusters_to_fasta(self, clusters_list: str | Path) -> None:
        """
        Extract sequences from cluster(s) and throw them in file(s)
        named cluster_X.faa in `output_dir/fastas/`.

        Parameter:
        -----------
        `clusters_list`: path to a file containing one cluster id (ex: cluster_0) per line.
        """
        self.output_dir.joinpath("fastas/").mkdir(parents=True, exist_ok=True)
        fasta_dir = self.output_dir.joinpath("fastas/")
        with open(clusters_list, "r", encoding="UTF-8") as clstr_list:
            for clstr in clstr_list:
                cluster_name = clstr.strip()
                filtered_df = self.seq_info_df.filter(pl.col("cluster") == cluster_name)
                output_file_name = fasta_dir.joinpath(f"{cluster_name}.faa")
                with open(output_file_name, "w", encoding="UTF-8") as fasta_output:
                    filtered_df.apply(lambda row: fasta_output.write(f"{row[0]}\n{row[1]}\n"))

    def _matrix_from_cdhit(self) -> None:
        """
        Parse .clstr file (output from `CD-HIT`).
        - count number of hits per sample in each cluster;
        - write the .tsv output in `output_dir/tsv/`;
        - adapted from original script by @Pier-Luc Plante.
        """

        self.output_dir.joinpath("tsv/").mkdir(parents=True, exist_ok=True)
        tsv_dir = self.output_dir.joinpath("tsv/")
        output_file = tsv_dir.joinpath(f"{self.clstr_base_name}.tsv")

        grouped_counts = self.seq_info_df.groupby(["cluster", "source"]).agg(pl.count().alias("count"))
        matrix = grouped_counts.pivot(values='count', index='cluster', columns='source').fill_null(0)
        matrix.write_csv(output_file, separator='\t')


    def featurecounts_matrix_from_cdhit(self, featurecounts_table: str | Path) -> None:
        """
        Count the number of features by clusters and write a matrix in the form of a .tsv
        file where each row is a cluster and each column a sample.

        Parameter:
        -----------
        `featurecounts_table`: path to a .tsv file containaing Id and abondance of multiple samples.
        """
        featurecounts_df = pl.read_csv(featurecounts_table, separator='\t')
        col_featurecounts_df = featurecounts_df.columns
        seq_clust_df = self.seq_info_df
        df_join = seq_clust_df.join(featurecounts_df, left_on="seq_id", right_on=col_featurecounts_df[0], how="left")
        cluster_counts = df_join.pivot(index="cluster", columns="source", values=col_featurecounts_df[1], aggregate_function="sum")

        output_file = self.output_dir.joinpath("clusters_abundance.tsv")
        cluster_counts.write_csv(output_file, separator='\t')


    def cluster_annotation(self, annotation_table: str | Path) -> None:
        """
        Annotate the cluster representant 

        Parameter:
        -----------
        `annotation_table`: path to a .tsv file containaing Id and abondance of multiple samples.
        """
        annotation_df = pl.read_csv(annotation_table, separator='\t')
        col_annotation_df = annotation_df.columns
        seq_clust_df = self.seq_info_df
        df_join = seq_clust_df.join(annotation_df, left_on="seq_id", right_on=col_annotation_df[0], how="left")
        df_join = df_join.filter(pl.col("representative"))
        df_join = df_join.select(pl.col("^(seq_id|cluster|(.*)_(name|id|E_value))$"))

        output_file = self.output_dir.joinpath("clusters_annotation.tsv")
        df_join.write_csv(output_file, separator='\t')


    def to_polars(self, tsv_table: str | Path) -> None:
        """
        Write a pickled polars table from tsv file.
        """
        table = pl.read_csv(tsv_table, separator='\t')
        with open(self.output_dir.joinpath(f"{tsv_table.split('/')[-1].split('.')[0]}.pkl"), "wb") as f_output:
            pickle.dump(table, f_output)


    def _clstr_rep_fasta(self) -> None:
        """
        Write all the clstr representatives to a fasta file.
        """
        raise NotImplementedError



if __name__ == '__main__':


    parser = ArgumentParser(description='''`clusterseq` is a command-line tool to extract
                                            information about a clustering experiment made
                                            with software `CD-HIT`.''')
    parser.add_argument('-c', '--clstr', help='''.clstr file (CD-HIT's output)''')
    parser.add_argument('-f', '--fasta', help='''source .faa file (CD-HIT's input)''')
    parser.add_argument('-o', '--output', help='''output directory''', default='./')
    parser.add_argument('-e', '--extract', help='''file containing clusters of interest
                                                   (1 per line) -> output: 1 fasta file
                                                   per cluster.''')
    parser.add_argument('-a', '--abundance', help='''featurecounts abundance table''')
    parser.add_argument('-n', '--annotation', help='''hmmer and diammond annotation in a tsv file
                                                    generated by hmmer_summary.py ''')
    parser.add_argument('-b', '--binAnnot', help='''bin annotation table''')
    parser.add_argument('-t', '--contigs2bins', help='''contigs2bins table''')

    args = parser.parse_args()


    exp = ClusterSeq(args.clstr, args.fasta, args.output)
    if args.extract:
        exp.clusters_to_fasta(args.extract)
    if args.abundance:
        exp.featurecounts_matrix_from_cdhit(args.abundance)
    if args.annotation:
        exp.cluster_annotation(args.annotation)
    if args.binAnnot:
        exp.to_polars(args.binAnnot)
    if args.contigs2bins:
        exp.to_polars(args.contigs2bins)
