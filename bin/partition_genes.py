#!/usr/bin/env python3
"""Split the proteins of one assembly between the final MAGs and everything else.

A protein belongs to a final MAG when the contig it was predicted on does: Prodigal and
Pyrodigal name each gene `<contig>_<n>`, so the contig is the gene id minus its last field.
"""

import argparse
import csv


def parse_args():
    parser = argparse.ArgumentParser(description="split proteins between final-MAG genes and the rest.")
    parser.add_argument("--faa", required=True, help="proteins predicted on the assembly.")
    parser.add_argument("--final-contig2mag", required=True, help="final_contig2mag.tsv written by summarize_mags.py.")
    parser.add_argument("--prefix", required=True, help="prefix of the two output files.")
    return parser.parse_args()


def contig_of(gene_id):
    return gene_id.rsplit("_", 1)[0]


def main():
    args = parse_args()

    with open(args.final_contig2mag) as handle:
        final_contigs = {row["contig"] for row in csv.DictReader(handle, delimiter="\t")}

    with open(args.faa) as faa, open(f"{args.prefix}.mag_genes.faa", "w") as mag_out, \
            open(f"{args.prefix}.other_genes.faa", "w") as other_out:
        out = other_out
        for line in faa:
            if line.startswith(">"):
                out = mag_out if contig_of(line[1:].split()[0]) in final_contigs else other_out
            out.write(line)


if __name__ == "__main__":
    main()
