#!/usr/bin/env python3
"""Give every input MAG its fate and collect the final, dereplicated MAG set.

    final        the representative of its dRep cluster (or the only MAG that passed)
    redundant    passed the quality filter, but another MAG of its cluster represents it
    low_quality  failed the quality filter, never dereplicated

dRep only runs when at least two MAGs pass the filter; without its tables, every MAG that
passed is final.
"""

import argparse
import csv
import os
import shutil
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="summarise MAG dereplication and collect the final MAGs.")
    parser.add_argument("--quality", required=True, help="mag_quality.tsv written by filter_mags.py.")
    parser.add_argument("--cdb", help="dRep data_tables/Cdb.csv (cluster of every genome).")
    parser.add_argument("--sdb", help="dRep data_tables/Sdb.csv (score of every genome).")
    parser.add_argument("--wdb", help="dRep data_tables/Wdb.csv (winner of every cluster).")
    parser.add_argument("--contig2mag", nargs="+", required=True, help="contig2mag.tsv tables written by prepare_assembly.py.")
    parser.add_argument("--mags", nargs="*", default=[], help="renamed MAG FASTA files.")
    parser.add_argument("--final-dir", default="final_mags", help="output directory for the final MAGs.")
    args = parser.parse_args()
    drep_tables = [args.cdb, args.sdb, args.wdb]
    if any(drep_tables) and not all(drep_tables):
        parser.error("--cdb, --sdb and --wdb go together")
    return args


def read_csv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle))


def genome_id(file_name):
    return file_name.removesuffix(".fa")


def main():
    args = parse_args()

    with open(args.quality) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    passing = [row["mag_id"] for row in rows if row["passes_filter"] == "yes"]

    cluster, score, representative = {}, {}, {}
    if args.cdb:
        cluster = {genome_id(r["genome"]): r["secondary_cluster"] for r in read_csv(args.cdb)}
        score = {genome_id(r["genome"]): r["score"] for r in read_csv(args.sdb)}
        winner_of = {r["cluster"]: genome_id(r["genome"]) for r in read_csv(args.wdb)}
        if set(cluster) != set(passing):
            sys.exit(f"error: dRep clustered {sorted(cluster)}, but the MAGs passing the filter are {sorted(passing)}")
        representative = {mag: winner_of[c] for mag, c in cluster.items()}
    elif len(passing) > 1:
        sys.exit(f"error: {len(passing)} MAGs passed the filter but no dRep table was given")
    else:
        representative = {mag: mag for mag in passing}

    columns = ["mag_id", "assembly", "original_bin", "n_contigs", "genome_size", "n50", "completeness",
               "contamination", "status", "secondary_cluster", "representative", "drep_score"]
    final = set()
    with open("mags_summary.tsv", "w") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            mag = row["mag_id"]
            if row["passes_filter"] != "yes":
                status = "low_quality"
            elif representative[mag] == mag:
                status = "final"
                final.add(mag)
            else:
                status = "redundant"
            row.update(status=status,
                       secondary_cluster=cluster.get(mag, ""),
                       representative=representative.get(mag, ""),
                       drep_score=score.get(mag, ""))
            out.write("\t".join(row[c] for c in columns) + "\n")

    contigs_of = {}
    for path in sorted(args.contig2mag):
        with open(path) as handle:
            for r in csv.DictReader(handle, delimiter="\t"):
                if r["mag_id"] in final:
                    contigs_of.setdefault(r["mag_id"], []).append(r["contig"])
    with open("final_contig2mag.tsv", "w") as out:
        out.write("contig\tmag_id\n")
        for mag in sorted(contigs_of):
            for contig in contigs_of[mag]:
                out.write(f"{contig}\t{mag}\n")

    # Copied rather than linked: a link to the absolute path of a staged input only resolves
    # where the task ran, which is not a given on cloud executors.
    os.makedirs(args.final_dir, exist_ok=True)
    files = {genome_id(os.path.basename(p)): p for p in args.mags}
    missing = sorted(final - set(files))
    if missing:
        sys.exit(f"error: no FASTA file for final MAGs: {', '.join(missing)}")
    for mag in sorted(final):
        shutil.copyfile(files[mag], os.path.join(args.final_dir, f"{mag}.fa"))


if __name__ == "__main__":
    main()
