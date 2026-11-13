#!/usr/bin/env python3
"""Group the dereplicated MAGs at a second, finer ANI level, from the comparisons dRep made.

dRep's secondary clustering is a hierarchical clustering of each primary cluster, cut at one ANI
threshold (the final MAGs, 95% by default). The same tree cut at a higher ANI gives finer groups
nested in dRep's clusters: 98% by default, the level at which dRep and inStrain consider genomes
distinct for read mapping.

The tree is rebuilt the way dRep builds it (drep.d_cluster.utils.make_linkage_Ndb): the ANI of a
comparison is set to 0 when its alignment coverage is at or below the coverage threshold, the two
directions are averaged, and the linkage runs on 1 - ANI over the genomes sorted by name. Cut at
dRep's own threshold, it must give back Cdb.csv exactly; the run stops otherwise, rather than
report groups from a tree dRep did not build.

    mag_levels.tsv  mag_id, secondary_cluster, fine_cluster, fine_representative

A fine group is named after its dRep cluster, `<secondary_cluster>.<n>`, numbered by the dRep
score of its representative, so that `.1` holds the final MAG.
"""

import argparse
import csv
import sys

from scipy.cluster.hierarchy import fcluster, linkage


def parse_args():
    parser = argparse.ArgumentParser(description="group dereplicated MAGs at a finer ANI level.")
    parser.add_argument("--ndb", required=True, help="dRep data_tables/Ndb.csv (secondary comparisons).")
    parser.add_argument("--cdb", required=True, help="dRep data_tables/Cdb.csv (cluster of every genome).")
    parser.add_argument("--sdb", required=True, help="dRep data_tables/Sdb.csv (score of every genome).")
    parser.add_argument("--fine-ani", type=float, required=True, help="ANI of the finer groups, e.g. 0.98.")
    parser.add_argument("--cov-thresh", type=float, default=0.1,
                        help="dRep's -nc, the coverage at or below which a comparison is ignored (dRep's default: 0.1).")
    return parser.parse_args()


def read_csv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle))


def genome_id(file_name):
    return file_name.removesuffix(".fa")


def partition(genomes, labels):
    """the groups of `genomes` that share a label, as a set of frozensets"""
    groups = {}
    for genome, label in zip(genomes, labels):
        groups.setdefault(label, set()).add(genome)
    return {frozenset(group) for group in groups.values()}


def cut_tree(genomes, ani, method, thresholds):
    """flat clusters of `genomes` at every distance threshold, from one average-linkage tree"""
    if len(genomes) == 1:
        return [[1] for _ in thresholds]
    condensed = []
    for i, a in enumerate(genomes):
        for b in genomes[i + 1:]:
            if (a, b) not in ani or (b, a) not in ani:
                sys.exit(f"error: Ndb.csv has no comparison of {a} with {b} in both directions")
            condensed.append(1 - (ani[(a, b)] + ani[(b, a)]) / 2)
    tree = linkage(condensed, method=method)
    return [list(fcluster(tree, threshold, criterion="distance")) for threshold in thresholds]


def main():
    args = parse_args()

    cdb = read_csv(args.cdb)
    if not cdb:
        sys.exit("error: Cdb.csv lists no genome")
    settings = {(row["threshold"], row["cluster_method"]) for row in cdb}
    if len(settings) != 1:
        sys.exit(f"error: Cdb.csv mixes clustering settings: {sorted(settings)}")
    threshold, method = settings.pop()
    threshold = float(threshold)
    fine_threshold = 1 - args.fine_ani
    if fine_threshold > threshold:
        sys.exit(f"error: --fine-ani {args.fine_ani} is below dRep's own ANI threshold ({1 - threshold:.4g}): "
                 "the finer level must be at least as strict")

    score = {row["genome"]: float(row["score"]) for row in read_csv(args.sdb)}
    unscored = sorted(row["genome"] for row in cdb if row["genome"] not in score)
    if unscored:
        sys.exit(f"error: Sdb.csv has no score for {', '.join(unscored)}")

    # (reference, querry) -> ANI, by primary cluster; the last value wins, as in dRep
    ani = {}
    for row in read_csv(args.ndb):
        value = 0.0 if float(row["alignment_coverage"]) <= args.cov_thresh else float(row["ani"])
        ani.setdefault(row["primary_cluster"], {})[(row["reference"], row["querry"])] = value

    members = {}
    for row in cdb:
        members.setdefault(row["primary_cluster"], []).append(row["genome"])
    secondary = {row["genome"]: row["secondary_cluster"] for row in cdb}

    rows = []
    for primary in sorted(members, key=int):
        genomes = sorted(members[primary])
        coarse, fine = cut_tree(genomes, ani.get(primary, {}), method, [threshold, fine_threshold])
        drep_groups = partition(genomes, [secondary[g] for g in genomes])
        if partition(genomes, coarse) != drep_groups:
            sys.exit(f"error: the tree rebuilt from Ndb.csv does not give back dRep's clusters of primary "
                     f"cluster {primary}; was dRep run with another -nc than --cov-thresh {args.cov_thresh}?")
        groups = partition(genomes, fine)
        if any(len({secondary[g] for g in group}) > 1 for group in groups):
            sys.exit(f"error: a fine group of primary cluster {primary} spans two dRep clusters")

        numbered = {}
        for group in sorted(groups, key=lambda g: (-max(score[m] for m in g), min(g))):
            cluster = secondary[next(iter(group))]
            numbered[cluster] = numbered.get(cluster, 0) + 1
            representative = min(group, key=lambda m: (-score[m], m))
            for genome in group:
                rows.append((genome_id(genome), cluster, f"{cluster}.{numbered[cluster]}", genome_id(representative)))

    with open("mag_levels.tsv", "w") as out:
        out.write("mag_id\tsecondary_cluster\tfine_cluster\tfine_representative\n")
        for row in sorted(rows):
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
