#!/usr/bin/env python3
"""Merge the three clustering steps into one non-redundant gene catalog.

    mag_cluster   final-MAG genes clustered among themselves
    mag_assigned  other genes absorbed by a final-MAG representative
    de_novo       the remaining genes clustered among themselves

The catalog holds one gene per cluster, its representative. The run stops here, before any
catalog file is published, if the clusters break any of these rules:

    - every gene id is unique across the assemblies
    - every gene is the member of exactly one cluster
    - every representative is a member of its own cluster
    - every final-MAG gene is represented by a final-MAG gene
"""

import argparse
import csv
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="build the gene catalog from the clustering steps.")
    parser.add_argument("--mag-clusters", help="mmseqs2 cluster table (representative, member) of the final-MAG genes.")
    parser.add_argument("--assigned", required=True, help="assigned.tsv written by assign_genes.py.")
    parser.add_argument("--other-clusters", help="mmseqs2 cluster table of the genes left unassigned.")
    parser.add_argument("--faa", nargs="+", required=True, help="proteins predicted on each assembly, <assembly>.faa.")
    parser.add_argument("--fna", nargs="+", required=True, help="gene nucleotide sequences, <assembly>.fna.")
    parser.add_argument("--gff", nargs="+", required=True, help="gene coordinates, <assembly>.gff.")
    parser.add_argument("--contig2mag", nargs="+", required=True, help="contig2mag.tsv tables written by prepare_assembly.py.")
    parser.add_argument("--mags-summary", required=True, help="mags_summary.tsv written by summarize_mags.py.")
    return parser.parse_args()


def contig_of(gene_id):
    return gene_id.rsplit("_", 1)[0]


def read_genes(faa_paths):
    """gene id -> assembly, in input order"""
    assembly_of = {}
    for path in sorted(faa_paths):
        assembly = os.path.basename(path).removesuffix(".faa")
        with open(path) as faa:
            for line in faa:
                if line.startswith(">"):
                    gene = line[1:].split()[0]
                    if gene in assembly_of:
                        sys.exit(f"error: gene id '{gene}' is predicted in both {assembly_of[gene]} and {assembly}")
                    assembly_of[gene] = assembly
    return assembly_of


def read_clusters(assembly_of, tables):
    rep_of, step_of = {}, {}
    for path, step in tables:
        if not path:
            continue
        with open(path) as handle:
            for line in handle:
                rep, member = line.rstrip("\n").split("\t")
                for gene in (rep, member):
                    if gene not in assembly_of:
                        sys.exit(f"error: {os.path.basename(path)} names gene '{gene}', which no assembly predicts")
                if member in rep_of:
                    sys.exit(f"error: gene '{member}' is a member of two clusters ({step_of[member]} and {step})")
                rep_of[member] = rep
                step_of[member] = step
    return rep_of, step_of


def write_records(paths, keep, output):
    """copy the FASTA records whose id is in `keep`; return how many were written"""
    written = 0
    with open(output, "w") as out:
        for path in sorted(paths):
            with open(path) as handle:
                copy = False
                for line in handle:
                    if line.startswith(">"):
                        copy = line[1:].split()[0] in keep
                        written += copy
                    if copy:
                        out.write(line)
    return written


def write_gff(paths, keep, output):
    written = 0
    with open(output, "w") as out:
        out.write("##gff-version 3\n")
        for path in sorted(paths):
            with open(path) as handle:
                for line in handle:
                    if line.startswith("#") or not line.strip():
                        continue
                    attributes = line.rstrip("\n").split("\t")[8]
                    gene = next((a[3:] for a in attributes.split(";") if a.startswith("ID=")), None)
                    if gene is None:
                        sys.exit(f"error: a feature of {os.path.basename(path)} has no ID attribute: {line.strip()}")
                    if gene in keep:
                        out.write(line)
                        written += 1
    return written


def main():
    args = parse_args()

    mag_of_contig = {}
    for path in args.contig2mag:
        with open(path) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                mag_of_contig[row["contig"]] = row["mag_id"]
    with open(args.mags_summary) as handle:
        status_of_mag = {row["mag_id"]: row["status"] for row in csv.DictReader(handle, delimiter="\t")}

    def mag_of(gene):
        return mag_of_contig.get(contig_of(gene), "")

    def status_of(gene):
        mag = mag_of(gene)
        return status_of_mag[mag] if mag else "unbinned"

    assembly_of = read_genes(args.faa)
    rep_of, step_of = read_clusters(assembly_of, [
        (args.mag_clusters, "mag_cluster"),
        (args.assigned, "mag_assigned"),
        (args.other_clusters, "de_novo"),
    ])

    unclustered = [gene for gene in assembly_of if gene not in rep_of]
    if unclustered:
        sys.exit(f"error: {len(unclustered)} gene(s) belong to no cluster, e.g. {unclustered[0]}")
    reps = set(rep_of.values())
    not_own = sorted(rep for rep in reps if rep_of[rep] != rep)
    if not_own:
        sys.exit(f"error: representative '{not_own[0]}' is itself a member of another cluster")
    unforced = [gene for gene in assembly_of if status_of(gene) == "final" and status_of(rep_of[gene]) != "final"]
    if unforced:
        gene = unforced[0]
        sys.exit(f"error: {len(unforced)} final-MAG gene(s) are represented by a gene outside the final MAGs, "
                 f"e.g. '{gene}' by '{rep_of[gene]}' ({status_of(rep_of[gene])})")

    with open("gene_clusters.tsv", "w") as out:
        out.write("representative\tmember\tmember_assembly\tmember_contig\tmember_mag\tmember_mag_status\tstep\n")
        for gene, assembly in assembly_of.items():
            out.write("\t".join([rep_of[gene], gene, assembly, contig_of(gene), mag_of(gene),
                                 status_of(gene), step_of[gene]]) + "\n")

    counts = {
        "faa": write_records(args.faa, reps, "gene_catalog.faa"),
        "fna": write_records(args.fna, reps, "gene_catalog.fna"),
        "gff": write_gff(args.gff, reps, "gene_catalog.gff"),
    }
    wrong = {kind: n for kind, n in counts.items() if n != len(reps)}
    if wrong:
        sys.exit(f"error: the catalog has {len(reps)} representatives, but wrote {wrong}")

    final_contigs = {contig for contig, mag in mag_of_contig.items() if status_of_mag[mag] == "final"}
    stats = [
        ("genes", len(assembly_of)),
        ("final_mag_genes", sum(status_of(g) == "final" for g in assembly_of)),
        ("representatives", len(reps)),
        ("final_mag_representatives", sum(status_of(r) == "final" for r in reps)),
        ("genes_absorbed_by_final_mag_representatives", sum(s == "mag_assigned" for s in step_of.values())),
        ("de_novo_clusters", sum(step_of[r] == "de_novo" for r in reps)),
        ("final_mag_contigs", len(final_contigs)),
        ("final_mag_contigs_in_catalog_gff", len({contig_of(r) for r in reps} & final_contigs)),
    ]
    with open("gene_catalog_stats.tsv", "w") as out:
        out.write("metric\tvalue\n")
        for metric, value in stats:
            out.write(f"{metric}\t{value}\n")


if __name__ == "__main__":
    main()
