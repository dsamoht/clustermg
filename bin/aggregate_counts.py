#!/usr/bin/env python3
"""Gather the counts of every read library, mapped once onto all the assemblies.

Every library was mapped once onto the contigs of every assembly, so each read sits on one copy
only: adding copies up counts it once. DNA is counted on contigs, RNA on genes.

Raw counts, one row per copy and library, only rows that are not 0:

    contig_dna_counts.tsv  contig, library: reads, aligned_bases, covered_bases
    link_dna_counts.tsv    link_id, library: the same, on the stretches link_contigs.py tied to a MAG
    mag_dna_profile.tsv    mag_id, library: CoverM on the MAG's own contigs, every MAG and library
    gene_rna_counts.tsv    gene, library: fragments, split between the genes a fragment overlaps

The main sums, for the users who do not query the database:

    final_mag_dna.tsv      final MAG x DNA library: its dRep cluster, every member's contigs and
                           linked stretches added up; mean_depth = aligned bases / its genome size
    fine_group_dna.tsv     the same for the finer groups, named after their representative
    gene_expression.tsv    representative x RNA library: its copies added up, and TPM on the
                           representative's length
    libraries.tsv          every library with its read totals: the reads it came with and, for RNA,
                           the rRNA reads SortMeRNA removed before mapping (both counting mates
                           apart), what was mapped (samtools flagstat) and what featureCounts
                           assigned (fragments)
"""

import argparse
import csv
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="gather the counts of every read library.")
    parser.add_argument("--libraries", required=True, help="library, sample, molecule, platform, layout, strandedness.")
    parser.add_argument("--flagstat", nargs="+", required=True, help="samtools flagstat -O tsv, <library>.flagstat.tsv.")
    parser.add_argument("--bedcov", nargs="*", default=[], help="samtools bedcov -c -j -d 1 of the DNA libraries, <library>.bedcov.tsv.")
    parser.add_argument("--featurecounts", nargs="*", default=[], help="featureCounts tables of the RNA libraries, <library>.featureCounts.tsv.")
    parser.add_argument("--coverm", help="coverm genome table of the DNA libraries.")
    parser.add_argument("--sortmerna-log", nargs="*", default=[], help="SortMeRNA logs of the RNA libraries, <library>.sortmerna.log.")
    parser.add_argument("--contig-links", required=True, help="contig_links.tsv written by link_contigs.py.")
    parser.add_argument("--contig2mag", nargs="+", required=True, help="contig2mag.tsv tables written by prepare_assembly.py.")
    parser.add_argument("--mags-summary", required=True, help="mags_summary.tsv written by summarize_mags.py.")
    parser.add_argument("--gene-clusters", required=True, help="gene_clusters.tsv written by build_catalog.py.")
    parser.add_argument("--saf", required=True, help="reference_genes.saf written by prepare_reference.py.")
    return parser.parse_args()


def fmt(value):
    return f"{value:.4f}".rstrip("0").rstrip(".")


def by_library(paths, suffix):
    return {os.path.basename(p).removesuffix(suffix): p for p in paths}


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_flagstat(path):
    values = {}
    with open(path) as handle:
        for line in handle:
            passed, _failed, what = line.rstrip("\n").split("\t")
            values[what] = passed
    return values.get("primary", ""), values.get("primary mapped", "")


def read_sortmerna_total(path):
    """the reads SortMeRNA was given, mates counted apart"""
    with open(path) as handle:
        for line in handle:
            if line.strip().startswith("Total reads ="):
                return int(line.split("=")[1])
    sys.exit(f"error: {path} has no 'Total reads' line")


def read_featurecounts_summary(path):
    with open(path) as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle][1:]
    return {status: float(value) for status, value in rows}


COVERM_METHODS = [(" Read Count", "reads"), (" Trimmed Mean", "trimmed_mean"),
                  (" Covered Fraction", "covered_fraction"), (" Mean", "mean")]


def read_coverm(path, libraries):
    """(mag, library) -> {reads, mean, trimmed_mean, covered_fraction}"""
    profile = {}
    with open(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        columns = []
        for column in header[1:]:
            suffix, method = next((s, m) for s, m in COVERM_METHODS if column.endswith(s))
            library = column.removesuffix(suffix)
            if library not in libraries:
                sys.exit(f"error: CoverM column '{column}' names no DNA library")
            columns.append((library, method))
        for row in reader:
            if row[0] == "unmapped":
                continue
            for (library, method), value in zip(columns, row[1:]):
                profile.setdefault((row[0], library), {})[method] = value
    return profile


def main():
    args = parse_args()

    libraries = {row["library"]: row for row in read_tsv(args.libraries)}
    dna = sorted(lib for lib, row in libraries.items() if row["molecule"] == "dna")
    rna = sorted(lib for lib, row in libraries.items() if row["molecule"] == "rna")
    flagstat = by_library(args.flagstat, ".flagstat.tsv")
    bedcov = by_library(args.bedcov, ".bedcov.tsv")
    featurecounts = by_library(args.featurecounts, ".featureCounts.tsv")
    sortmerna = by_library(args.sortmerna_log, ".sortmerna.log")
    if sortmerna and sorted(sortmerna) != rna:
        sys.exit(f"error: the SortMeRNA logs are for libraries {sorted(sortmerna)}, expected {rna}")
    for kind, files, expected in [("flagstat", flagstat, sorted(libraries)), ("bedcov", bedcov, dna),
                                  ("featureCounts", featurecounts, rna)]:
        if sorted(files) != expected:
            sys.exit(f"error: the {kind} tables are for libraries {sorted(files)}, expected {expected}")

    mags = {row["mag_id"]: row for row in read_tsv(args.mags_summary)}
    passing = {mag: row for mag, row in mags.items() if row["status"] in ("final", "redundant")}
    mag_of = {}
    for path in args.contig2mag:
        for row in read_tsv(path):
            mag_of[row["contig"]] = row["mag_id"]
    link_mag = {row["link_id"]: row["mag_id"] for row in read_tsv(args.contig_links)}

    # DNA: contigs and linked stretches, then every passing MAG's own share
    member_counts = {}  # (passing mag, library) -> [reads, aligned bases]
    with open("contig_dna_counts.tsv", "w") as contigs_out, open("link_dna_counts.tsv", "w") as links_out:
        contigs_out.write("contig\tlibrary\treads\taligned_bases\tcovered_bases\n")
        links_out.write("link_id\tlibrary\treads\taligned_bases\tcovered_bases\n")
        for library in dna:
            with open(bedcov[library]) as handle:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    contig, name = fields[0], fields[3]
                    bases, covered, reads = int(fields[4]), int(fields[5]), int(fields[6])
                    if name == "contig":
                        mag = mag_of.get(contig) if mag_of.get(contig) in passing else None
                        if reads or bases:
                            contigs_out.write(f"{contig}\t{library}\t{reads}\t{bases}\t{covered}\n")
                    else:
                        mag = link_mag.get(name)
                        if mag is None:
                            sys.exit(f"error: region '{name}' of {library} is in no link of contig_links.tsv")
                        if reads or bases:
                            links_out.write(f"{name}\t{library}\t{reads}\t{bases}\t{covered}\n")
                    if mag is not None:
                        total = member_counts.setdefault((mag, library), [0, 0])
                        total[0] += reads
                        total[1] += bases

    for level, output in [("representative", "final_mag_dna.tsv"), ("fine_representative", "fine_group_dna.tsv")]:
        groups = sorted({row[level] for row in passing.values()})
        totals = {}
        for mag, row in passing.items():
            for library in dna:
                reads, bases = member_counts.get((mag, library), [0, 0])
                total = totals.setdefault((row[level], library), [0, 0])
                total[0] += reads
                total[1] += bases
        with open(output, "w") as out:
            out.write("mag_id\tlibrary\treads\taligned_bases\tmean_depth\n")
            for group in groups:
                size = int(mags[group]["genome_size"])
                for library in dna:
                    reads, bases = totals[(group, library)]
                    out.write(f"{group}\t{library}\t{reads}\t{bases}\t{fmt(bases / size)}\n")

    profile = read_coverm(args.coverm, set(dna)) if args.coverm else {}
    with open("mag_dna_profile.tsv", "w") as out:
        out.write("mag_id\tlibrary\treads\tmean\ttrimmed_mean\tcovered_fraction\n")
        for (mag, library), values in sorted(profile.items()):
            if mag not in mags:
                sys.exit(f"error: CoverM profiled '{mag}', which is no MAG of mags_summary.tsv")
            out.write("\t".join([mag, library, values["reads"], values["mean"], values["trimmed_mean"],
                                 values["covered_fraction"]]) + "\n")

    # RNA: every gene copy, then its cluster's representative
    representative = {row["member"]: row["representative"] for row in read_tsv(args.gene_clusters)}
    length = {row["GeneID"]: int(row["End"]) - int(row["Start"]) + 1 for row in read_tsv(args.saf)}
    assigned = {}
    expression = {}  # library -> representative -> fragments
    with open("gene_rna_counts.tsv", "w") as out:
        out.write("gene\tlibrary\treads\n")
        for library in rna:
            summary = read_featurecounts_summary(featurecounts[library] + ".summary")
            assigned[library] = summary
            total, fractional = 0.0, 0
            sums = expression.setdefault(library, {})
            with open(featurecounts[library]) as handle:
                lines = (line for line in handle if not line.startswith("#"))
                next(lines)
                for line in lines:
                    fields = line.rstrip("\n").split("\t")
                    count = float(fields[-1])
                    if not count:
                        continue
                    gene = fields[0]
                    if gene not in representative:
                        sys.exit(f"error: gene '{gene}' counted for {library} is in no cluster of gene_clusters.tsv")
                    out.write(f"{gene}\t{library}\t{fmt(count)}\n")
                    total += count
                    fractional += not count.is_integer()
                    rep = representative[gene]
                    sums[rep] = sums.get(rep, 0.0) + count
            # featureCounts writes a split fragment's shares with two decimals: each count off by
            # up to 0.005
            if abs(total - summary["Assigned"]) > 0.005 * fractional + 1e-6 * summary["Assigned"] + 1e-9:
                sys.exit(f"error: the gene counts of {library} add up to {total}, but featureCounts assigned "
                         f"{summary['Assigned']:.0f} fragments: a fragment would be counted twice or lost")

    with open("gene_expression.tsv", "w") as out:
        out.write("representative\tlibrary\treads\ttpm\n")
        for library in rna:
            sums = expression[library]
            per_kb = {rep: count / (length[rep] / 1000) for rep, count in sums.items()}
            scale = sum(per_kb.values())
            for rep in sorted(sums):
                out.write(f"{rep}\t{library}\t{fmt(sums[rep])}\t{fmt(per_kb[rep] / scale * 1e6)}\n")

    with open("libraries.tsv", "w") as out:
        out.write("library\tsample\tmolecule\tplatform\tlayout\tstrandedness\tinput_reads\trrna_reads\t"
                  "primary_reads\tprimary_mapped\tassigned\tunassigned_no_features\n")
        for library in sorted(libraries):
            row = libraries[library]
            reads, mapped = read_flagstat(flagstat[library])
            # every read that reached the mapper is in the BAM, mapped or not: what SortMeRNA was
            # given and did not pass on is rRNA
            if library in sortmerna:
                total = read_sortmerna_total(sortmerna[library])
                input_reads, rrna = str(total), str(total - int(reads))
            else:
                input_reads, rrna = reads, ""
            summary = assigned.get(library, {})
            counted = [fmt(summary[k]) if k in summary else "" for k in ("Assigned", "Unassigned_NoFeatures")]
            out.write("\t".join([library, row["sample"], row["molecule"], row["platform"], row["layout"],
                                 row["strandedness"], input_reads, rrna, reads, mapped, *counted]) + "\n")


if __name__ == "__main__":
    main()
