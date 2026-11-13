#!/usr/bin/env python3
"""Tie the stretches of unbinned and low-quality contigs to the MAG they are a copy of.

A genome assembled in several assemblies is binned in some and left unbinned, or in a bin that
failed the filter, in others. Every read library is mapped onto all the assemblies, so some of
its reads land on those copies. A stretch of such a contig that aligns to a passing MAG (final
or redundant) at `--min-identity` over `--min-length` of the contig is counted as part of that
MAG. Identity is the matching bases over the alignment length, gaps included (PAF columns 10 and
11). A link only says that the stretch is that similar: a stretch of a related strain, or of a
sequence several genomes share (a plasmid, a transposon, an rRNA operon), links as well.

Where a stretch aligns to several MAGs, minimap2 (`--secondary=no`) reports only its
best-scoring alignment, so the choice between MAGs is mostly minimap2's: its score weighs the
matches against the mismatches and gaps. The alignments it does report can still overlap on the
contig, by less than half of the shorter one (minimap2's mask level). Those are settled here:
the alignment of higher identity, then the longer, is kept, and the others are dropped whole.

The contigs of passing MAGs are never linked: they belong to their own MAG.

    contig_links.tsv  link_id, contig, start, end, mag_id, identity, aligned_length
                      (start 0-based, end exclusive, as in BED; link_id is contig:start+1-end)
    regions.bed       contigs.bed, then one row per link: the regions DNA is counted on
"""

import argparse
import csv
import shutil


def parse_args():
    parser = argparse.ArgumentParser(description="link stretches of unbinned contigs to the MAG they copy.")
    parser.add_argument("--paf", help="minimap2 alignments of the other contigs (query) to the passing-MAG contigs (target).")
    parser.add_argument("--contigs-bed", required=True, help="contigs.bed written by prepare_reference.py.")
    parser.add_argument("--contig2mag", nargs="+", required=True, help="contig2mag.tsv tables written by prepare_assembly.py.")
    parser.add_argument("--min-identity", type=float, required=True, help="minimum alignment identity, e.g. 0.98.")
    parser.add_argument("--min-length", type=int, required=True, help="minimum aligned length on the contig (bp).")
    return parser.parse_args()


def candidates(paf_path, mag_of, min_identity, min_length):
    """contig -> [(identity, length, mag, start, end)] of the alignments above the thresholds"""
    hits = {}
    with open(paf_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            contig, start, end, target = fields[0], int(fields[2]), int(fields[3]), fields[5]
            identity = int(fields[9]) / int(fields[10])
            if identity < min_identity or end - start < min_length:
                continue
            hits.setdefault(contig, []).append((identity, end - start, mag_of[target], start, end))
    return hits


def main():
    args = parse_args()

    mag_of = {}
    for path in args.contig2mag:
        with open(path) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                mag_of[row["contig"]] = row["mag_id"]

    links = []
    if args.paf:
        for contig, hits in candidates(args.paf, mag_of, args.min_identity, args.min_length).items():
            taken = []
            for identity, length, mag, start, end in sorted(hits, key=lambda h: (-h[0], -h[1], h[2], h[3])):
                if any(start < t_end and t_start < end for t_start, t_end in taken):
                    continue
                taken.append((start, end))
                links.append((contig, start, end, mag, identity, length))

    with open("contig_links.tsv", "w") as out:
        out.write("link_id\tcontig\tstart\tend\tmag_id\tidentity\taligned_length\n")
        for contig, start, end, mag, identity, length in sorted(links):
            out.write(f"{contig}:{start + 1}-{end}\t{contig}\t{start}\t{end}\t{mag}\t{identity:.4f}\t{length}\n")

    shutil.copyfile(args.contigs_bed, "regions.bed")
    with open("regions.bed", "a") as out:
        for contig, start, end, *_ in sorted(links):
            out.write(f"{contig}\t{start}\t{end}\t{contig}:{start + 1}-{end}\n")


if __name__ == "__main__":
    main()
