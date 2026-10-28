#!/usr/bin/env python3
"""Keep the MAGs that are good enough to be dereplicated, and write dRep's inputs for them.

The thresholds are the ones dRep is given as well, so dRep's own filter has nothing left to
remove; applying them here first is what tells the pipeline whether dRep has anything to do.
"""

import argparse
import csv
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="filter MAGs on completeness, contamination and length.")
    parser.add_argument("--min-completeness", type=float, required=True, help="minimum CheckM2 completeness (%%).")
    parser.add_argument("--max-contamination", type=float, required=True, help="maximum CheckM2 contamination (%%).")
    parser.add_argument("--min-length", type=int, required=True, help="minimum genome size (bp).")
    parser.add_argument("mag_info", nargs="+", help="mag_info.tsv tables written by prepare_assembly.py.")
    return parser.parse_args()


def passes(row, args):
    try:
        completeness = float(row["completeness"])
        contamination = float(row["contamination"])
    except ValueError:
        # a MAG CheckM2 could not assess is not dereplicated
        return False
    return (completeness >= args.min_completeness
            and contamination <= args.max_contamination
            and int(row["genome_size"]) >= args.min_length)


def main():
    args = parse_args()

    rows, columns, origin = [], None, {}
    for path in sorted(args.mag_info):
        with open(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = columns or reader.fieldnames
            for row in reader:
                # MAG ids name the files dRep sees, so two assemblies must not produce the same one
                if row["mag_id"] in origin:
                    sys.exit(f"error: MAG id '{row['mag_id']}' comes from assemblies {origin[row['mag_id']]} and {row['assembly']}")
                origin[row["mag_id"]] = row["assembly"]
                rows.append(row)
    rows.sort(key=lambda row: row["mag_id"])

    with open("mag_quality.tsv", "w") as quality, open("genome_info.csv", "w") as info, \
            open("passing_mags.txt", "w") as passing:
        header = columns + ["passes_filter"]
        quality.write("\t".join(header) + "\n")
        info.write("genome,completeness,contamination\n")
        for row in rows:
            ok = passes(row, args)
            quality.write("\t".join([row[c] for c in header[:-1]] + ["yes" if ok else "no"]) + "\n")
            if ok:
                info.write(f"{row['mag_id']}.fa,{row['completeness']},{row['contamination']}\n")
                passing.write(f"{row['mag_id']}.fa\n")


if __name__ == "__main__":
    main()
