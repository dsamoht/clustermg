#!/usr/bin/env python3
"""Rename the contigs and MAGs of one assembly and check its bins against it.

Contig names collide across assemblies (MEGAHIT's `k141_1`, Flye's `contig_1`), so every
contig becomes `<name>_<first word of its header>`. A MAG becomes `<name>_<bin file name>`,
unless the file name already starts with `<name>_`, as mag-ont's `<group>_bin_001.fa` does.

The MAG contigs must be the assembly's contigs, sequence for sequence: genes are predicted
on the assembly, and a final MAG's genes are found through the contigs it is made of.
"""

import argparse
import csv
import gzip
import hashlib
import os
import sys

FASTA_EXTENSIONS = (".fa", ".fasta", ".fna")
# IUPAC nucleotide codes plus gaps: anything else, protein letters included, is refused
DNA_ALPHABET = "ACGTUNRYKMSWBDHV-"
NOT_DNA = str.maketrans("", "", DNA_ALPHABET)


def parse_args():
    parser = argparse.ArgumentParser(description="rename the contigs and MAGs of an assembly and check its bins.")
    parser.add_argument("--name", required=True, help="assembly name, the prefix of every contig and MAG.")
    parser.add_argument("--assembly", required=True, help="assembly FASTA, optionally gzipped.")
    parser.add_argument("--bins", help="directory of bin FASTA files, optionally gzipped.")
    parser.add_argument("--checkm2", help="CheckM2 quality_report.tsv covering the bins.")
    parser.add_argument("--mag-dir", default="mags", help="output directory for the renamed MAGs.")
    args = parser.parse_args()
    if bool(args.bins) != bool(args.checkm2):
        parser.error("--bins and --checkm2 go together")
    return args


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_fasta(path):
    """yield (contig id, full header, sequence) from a fasta file, gzipped or not."""
    header, seq = None, []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header.split()[0], header, "".join(seq)
                header, seq = line[1:], []
            else:
                seq.append(line)
    if header is not None:
        yield header.split()[0], header, "".join(seq)


def write_record(handle, header, seq, width=60):
    handle.write(f">{header}\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i:i + width] + "\n")


def digest(seq):
    return hashlib.md5(seq.upper().encode()).digest()


def bin_basename(file_name):
    stem = file_name[:-3] if file_name.endswith(".gz") else file_name
    return os.path.splitext(stem)[0]


def mag_id_of(name, basename):
    return basename if basename.startswith(f"{name}_") else f"{name}_{basename}"


def n50(lengths):
    total, running = sum(lengths), 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running * 2 >= total:
            return length
    return 0


def list_bins(bins_dir):
    files = sorted(
        f for f in os.listdir(bins_dir)
        if f.lower().removesuffix(".gz").endswith(FASTA_EXTENSIONS)
    )
    if not files:
        sys.exit(f"error: no FASTA file ({', '.join(FASTA_EXTENSIONS)}, optionally .gz) in bins directory {bins_dir}")
    return [os.path.join(bins_dir, f) for f in files]


def read_bins(name, bin_paths):
    """map every binned contig to its MAG and sequence digest, and describe every MAG."""
    contig_to_bin = {}
    mags = []
    seen_ids = {}
    for path in bin_paths:
        basename = bin_basename(os.path.basename(path))
        mag_id = mag_id_of(name, basename)
        if mag_id in seen_ids:
            sys.exit(f"error: bins {seen_ids[mag_id]} and {os.path.basename(path)} would both become MAG '{mag_id}'")
        seen_ids[mag_id] = os.path.basename(path)

        lengths = []
        for contig, _header, seq in parse_fasta(path):
            if contig in contig_to_bin:
                other = contig_to_bin[contig][0]
                if other == mag_id:
                    sys.exit(f"error: contig '{contig}' is listed twice in bin {os.path.basename(path)}")
                sys.exit(f"error: contig '{contig}' is in two bins, {seen_ids[other]} and {os.path.basename(path)}")
            contig_to_bin[contig] = (mag_id, digest(seq))
            lengths.append(len(seq))
        if not lengths:
            sys.exit(f"error: bin {os.path.basename(path)} has no sequence")
        mags.append({
            "mag_id": mag_id,
            "basename": basename,
            "path": path,
            "n_contigs": len(lengths),
            "genome_size": sum(lengths),
            "n50": n50(lengths),
        })
    return contig_to_bin, mags


def read_checkm2(path, mags):
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing_cols = {"Name", "Completeness", "Contamination"} - set(reader.fieldnames or [])
        if missing_cols:
            sys.exit(f"error: {path} is not a CheckM2 quality report, it lacks {', '.join(sorted(missing_cols))}")
        quality = {row["Name"]: row for row in reader}
    missing = [mag["basename"] for mag in mags if mag["basename"] not in quality]
    if missing:
        sys.exit(f"error: no CheckM2 row in {path} for bins: {', '.join(missing)}")
    for mag in mags:
        row = quality[mag["basename"]]
        mag["completeness"] = row["Completeness"]
        mag["contamination"] = row["Contamination"]


def rename_assembly(name, assembly, contig_to_bin):
    seen = set()
    matched = set()
    with open(f"{name}.contigs.fa", "w") as contigs_out, open(f"{name}.contig_names.tsv", "w") as names_out:
        names_out.write("contig\toriginal_contig\n")
        for contig, _header, seq in parse_fasta(assembly):
            if contig in seen:
                sys.exit(f"error: contig name '{contig}' appears twice in assembly {assembly}")
            seen.add(contig)
            if seq.upper().translate(NOT_DNA):
                sys.exit(f"error: assembly {assembly} is not DNA: contig '{contig}' holds non-nucleotide characters")
            if contig in contig_to_bin:
                if contig_to_bin[contig][1] != digest(seq):
                    sys.exit(f"error: contig '{contig}' of MAG {contig_to_bin[contig][0]} differs from the same contig in assembly {assembly}")
                matched.add(contig)
            write_record(contigs_out, f"{name}_{contig}", seq)
            names_out.write(f"{name}_{contig}\t{contig}\n")
    missing = [contig for contig in contig_to_bin if contig not in matched]
    if missing:
        shown = ", ".join(missing[:5]) + (f" and {len(missing) - 5} more" if len(missing) > 5 else "")
        sys.exit(f"error: {len(missing)} bin contig(s) missing from assembly {assembly}: {shown}")


def write_mags(name, mags, mag_dir):
    os.makedirs(mag_dir, exist_ok=True)
    with open(f"{name}.contig2mag.tsv", "w") as c2m:
        c2m.write("contig\tmag_id\n")
        for mag in mags:
            with open(os.path.join(mag_dir, f"{mag['mag_id']}.fa"), "w") as out:
                for contig, _header, seq in parse_fasta(mag["path"]):
                    write_record(out, f"{name}_{contig}", seq)
                    c2m.write(f"{name}_{contig}\t{mag['mag_id']}\n")


def write_mag_info(name, mags):
    columns = ["mag_id", "assembly", "original_bin", "n_contigs", "genome_size", "n50", "completeness", "contamination"]
    with open(f"{name}.mag_info.tsv", "w") as out:
        out.write("\t".join(columns) + "\n")
        for mag in mags:
            row = [mag["mag_id"], name, os.path.basename(mag["path"]), mag["n_contigs"], mag["genome_size"],
                   mag["n50"], mag["completeness"], mag["contamination"]]
            out.write("\t".join(map(str, row)) + "\n")


def main():
    args = parse_args()

    contig_to_bin, mags = {}, []
    if args.bins:
        contig_to_bin, mags = read_bins(args.name, list_bins(args.bins))
        read_checkm2(args.checkm2, mags)

    rename_assembly(args.name, args.assembly, contig_to_bin)
    write_mags(args.name, mags, args.mag_dir)
    write_mag_info(args.name, mags)


if __name__ == "__main__":
    main()
