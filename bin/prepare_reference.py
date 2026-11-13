#!/usr/bin/env python3
"""Write the reference every read library is mapped to: the contigs of all the assemblies.

Every library is mapped once onto all the assemblies together, so that a read lands on one
copy only, the closest; adding the copies of a genome or a gene up then counts it once.

    reference.fa         every renamed contig of every assembly
    contigs.bed          every contig, whole, named `contig`: the regions DNA is counted on
    link_target.fa       the contigs of the MAGs that passed the filter (final and redundant)
    link_query.fa        every other contig: unbinned, or of a low-quality MAG
    reference_genes.saf  every predicted gene, for featureCounts

The stretches of the query contigs that align to a passing MAG are counted as part of it
(link_contigs.py): a genome's reads are not lost to a copy of it that no bin holds.
"""

import argparse
import csv
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="write the mapping reference of every read library.")
    parser.add_argument("--contigs", nargs="+", required=True, help="renamed contigs, <name>.contigs.fa.")
    parser.add_argument("--gff", nargs="+", required=True, help="predicted genes of every assembly, <name>.gff.")
    parser.add_argument("--contig2mag", nargs="+", required=True, help="contig2mag.tsv tables written by prepare_assembly.py.")
    parser.add_argument("--mags-summary", required=True, help="mags_summary.tsv written by summarize_mags.py.")
    return parser.parse_args()


def passing_contigs(contig2mag_paths, summary_path):
    with open(summary_path) as handle:
        passing = {row["mag_id"] for row in csv.DictReader(handle, delimiter="\t")
                   if row["status"] in ("final", "redundant")}
    contigs = set()
    for path in contig2mag_paths:
        with open(path) as handle:
            contigs.update(row["contig"] for row in csv.DictReader(handle, delimiter="\t") if row["mag_id"] in passing)
    return contigs


def write_reference(contig_paths, target_contigs):
    """copy every contig to the reference and to its side of the link alignment; return the contig ids"""
    seen = set()
    with open("reference.fa", "w") as reference, open("link_target.fa", "w") as target, \
            open("link_query.fa", "w") as query, open("contigs.bed", "w") as bed:
        for path in sorted(contig_paths):
            with open(path) as handle:
                contig, length, side = None, 0, None
                for line in handle:
                    if line.startswith(">"):
                        if contig is not None:
                            bed.write(f"{contig}\t0\t{length}\tcontig\n")
                        contig = line[1:].split()[0]
                        if contig in seen:
                            sys.exit(f"error: contig '{contig}' is in two assemblies")
                        seen.add(contig)
                        length = 0
                        side = target if contig in target_contigs else query
                    else:
                        length += len(line.strip())
                    reference.write(line)
                    side.write(line)
                if contig is not None:
                    bed.write(f"{contig}\t0\t{length}\tcontig\n")
    return seen


def write_saf(gff_paths, contigs):
    with open("reference_genes.saf", "w") as out:
        out.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for path in sorted(gff_paths):
            with open(path) as handle:
                for line in handle:
                    if line.startswith("#") or not line.strip():
                        continue
                    fields = line.rstrip("\n").split("\t")
                    gene = next((a[3:] for a in fields[8].split(";") if a.startswith("ID=")), None)
                    if gene is None:
                        sys.exit(f"error: a feature of {os.path.basename(path)} has no ID attribute: {line.strip()}")
                    if fields[0] not in contigs:
                        sys.exit(f"error: gene '{gene}' of {os.path.basename(path)} is on contig '{fields[0]}', which no assembly holds")
                    out.write("\t".join([gene, fields[0], fields[3], fields[4], fields[6]]) + "\n")


def main():
    args = parse_args()
    contigs = write_reference(args.contigs, passing_contigs(args.contig2mag, args.mags_summary))
    write_saf(args.gff, contigs)


if __name__ == "__main__":
    main()
