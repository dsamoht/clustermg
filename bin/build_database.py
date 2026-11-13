#!/usr/bin/env python3
"""Load every result of the run into one DuckDB database.

Tables hold the facts, once each; anything that can be derived is a view, so nothing is stored
twice. Foreign keys between tables are enforced as rows are loaded. The links inside a table
(MAG -> representative, gene -> representative) cannot be: DuckDB checks a foreign key against
the rows already in the table, and a final MAG is its own representative. They are checked
after loading instead, with the counts of every table, in database_checks.tsv, and the run
stops on any row that points nowhere.

    clustermg.duckdb    the database
    parquet/            every table as Parquet, with schema.sql and load.sql to rebuild it
    database_checks.tsv rows per table, and broken links (all 0)

Without read libraries, the count tables are created empty.
"""

import argparse
import csv
import os
import sys

import duckdb

SCHEMA = """
CREATE TABLE assembly (
    assembly_id VARCHAR PRIMARY KEY, kind VARCHAR, source VARCHAR);
CREATE TABLE sample (
    sample_id VARCHAR PRIMARY KEY);
CREATE TABLE assembly_sample (
    assembly_id VARCHAR REFERENCES assembly (assembly_id), sample_id VARCHAR REFERENCES sample (sample_id),
    PRIMARY KEY (assembly_id, sample_id));
CREATE TABLE library (
    library_id VARCHAR PRIMARY KEY, sample_id VARCHAR REFERENCES sample (sample_id), molecule VARCHAR,
    platform VARCHAR, layout VARCHAR, strandedness VARCHAR, input_reads BIGINT, rrna_reads BIGINT,
    primary_reads BIGINT, primary_mapped BIGINT, assigned DOUBLE, unassigned_no_features DOUBLE);
CREATE TABLE mag (
    mag_id VARCHAR PRIMARY KEY, assembly_id VARCHAR REFERENCES assembly (assembly_id), original_bin VARCHAR,
    n_contigs INTEGER, genome_size BIGINT, n50 BIGINT, completeness DOUBLE, contamination DOUBLE,
    status VARCHAR, secondary_cluster VARCHAR, representative VARCHAR, drep_score DOUBLE,
    fine_cluster VARCHAR, fine_representative VARCHAR);
CREATE TABLE contig (
    contig_id VARCHAR PRIMARY KEY, assembly_id VARCHAR REFERENCES assembly (assembly_id),
    original_contig VARCHAR, length BIGINT, mag_id VARCHAR REFERENCES mag (mag_id));
CREATE TABLE gene (
    gene_id VARCHAR PRIMARY KEY, contig_id VARCHAR REFERENCES contig (contig_id), start BIGINT, stop BIGINT,
    strand VARCHAR, partial VARCHAR, representative VARCHAR, step VARCHAR);
CREATE TABLE catalog_sequence (
    gene_id VARCHAR PRIMARY KEY REFERENCES gene (gene_id), protein VARCHAR, nucleotide VARCHAR);
CREATE TABLE contig_link (
    link_id VARCHAR PRIMARY KEY, contig_id VARCHAR REFERENCES contig (contig_id), start BIGINT, stop BIGINT,
    mag_id VARCHAR REFERENCES mag (mag_id), identity DOUBLE, aligned_length BIGINT);
CREATE TABLE contig_dna_counts (
    contig_id VARCHAR REFERENCES contig (contig_id), library_id VARCHAR REFERENCES library (library_id),
    reads BIGINT, aligned_bases BIGINT, covered_bases BIGINT, PRIMARY KEY (contig_id, library_id));
CREATE TABLE link_dna_counts (
    link_id VARCHAR REFERENCES contig_link (link_id), library_id VARCHAR REFERENCES library (library_id),
    reads BIGINT, aligned_bases BIGINT, covered_bases BIGINT, PRIMARY KEY (link_id, library_id));
CREATE TABLE mag_dna_profile (
    mag_id VARCHAR REFERENCES mag (mag_id), library_id VARCHAR REFERENCES library (library_id),
    reads BIGINT, mean DOUBLE, trimmed_mean DOUBLE, covered_fraction DOUBLE, PRIMARY KEY (mag_id, library_id));
CREATE TABLE gene_rna_counts (
    gene_id VARCHAR REFERENCES gene (gene_id), library_id VARCHAR REFERENCES library (library_id),
    reads DOUBLE, PRIMARY KEY (gene_id, library_id));
"""

VIEWS = """
-- The dereplicated MAGs, and the genes of the catalog
CREATE VIEW final_mag AS SELECT * FROM mag WHERE status = 'final';
CREATE VIEW catalog_gene AS
    SELECT g.*, s.protein, s.nucleotide FROM gene g JOIN catalog_sequence s USING (gene_id);

-- Where every gene sits: contig, assembly, MAG, finer group and final MAG
CREATE VIEW gene_context AS
    SELECT g.gene_id, g.representative, g.step, g.contig_id, c.assembly_id, a.kind AS assembly_kind,
           c.mag_id, m.status AS mag_status, m.fine_representative AS fine_group, m.representative AS final_mag
    FROM gene g JOIN contig c USING (contig_id) JOIN assembly a USING (assembly_id)
    LEFT JOIN mag m ON m.mag_id = c.mag_id;

-- DNA of every passing MAG: its own contigs and the stretches of other contigs linked to it
CREATE VIEW mag_dna_member AS
    SELECT c.mag_id, d.library_id, d.reads, d.aligned_bases
    FROM contig_dna_counts d JOIN contig c USING (contig_id) JOIN mag m ON m.mag_id = c.mag_id
    WHERE m.status IN ('final', 'redundant')
    UNION ALL
    SELECT l.mag_id, d.library_id, d.reads, d.aligned_bases
    FROM link_dna_counts d JOIN contig_link l USING (link_id);

-- A final MAG holds its whole dRep cluster; a finer group, its members. Each read is on one copy
-- only, so the members add up; mean depth is over the representative's genome size.
CREATE VIEW final_mag_dna AS
    SELECT m.representative AS mag_id, x.library_id, SUM(x.reads) AS reads,
           SUM(x.aligned_bases) AS aligned_bases, SUM(x.aligned_bases) / ANY_VALUE(r.genome_size) AS mean_depth
    FROM mag_dna_member x JOIN mag m USING (mag_id) JOIN mag r ON r.mag_id = m.representative
    GROUP BY m.representative, x.library_id;
CREATE VIEW fine_group_dna AS
    SELECT m.fine_representative AS mag_id, x.library_id, SUM(x.reads) AS reads,
           SUM(x.aligned_bases) AS aligned_bases, SUM(x.aligned_bases) / ANY_VALUE(r.genome_size) AS mean_depth
    FROM mag_dna_member x JOIN mag m USING (mag_id) JOIN mag r ON r.mag_id = m.fine_representative
    GROUP BY m.fine_representative, x.library_id;

-- RNA of every catalog gene: its copies added up, and TPM on the representative's length
CREATE VIEW gene_expression AS
    WITH per_gene AS (
        SELECT g.representative AS gene_id, c.library_id, SUM(c.reads) AS reads
        FROM gene_rna_counts c JOIN gene g USING (gene_id) GROUP BY g.representative, c.library_id),
    per_kb AS (
        SELECT p.*, p.reads / ((r.stop - r.start + 1) / 1000.0) AS rpk
        FROM per_gene p JOIN gene r ON r.gene_id = p.gene_id)
    SELECT gene_id, library_id, reads, rpk / SUM(rpk) OVER (PARTITION BY library_id) * 1e6 AS tpm FROM per_kb;

-- The passing MAG every gene copy belongs to: by its contig, or by a linked stretch it sits in
CREATE VIEW gene_mag AS
    SELECT g.gene_id, c.mag_id
    FROM gene g JOIN contig c USING (contig_id) JOIN mag m ON m.mag_id = c.mag_id
    WHERE m.status IN ('final', 'redundant')
    UNION ALL
    SELECT g.gene_id, l.mag_id
    FROM gene g JOIN contig_link l ON l.contig_id = g.contig_id AND (g.start + g.stop) / 2 BETWEEN l.start AND l.stop;
CREATE VIEW mag_expression AS
    SELECT m.representative AS mag_id, c.library_id, SUM(c.reads) AS reads
    FROM gene_rna_counts c JOIN gene_mag gm USING (gene_id) JOIN mag m ON m.mag_id = gm.mag_id
    GROUP BY m.representative, c.library_id;

-- DNA depth of every catalog gene: each copy takes the depth of its contig, the copies add up
CREATE VIEW gene_dna_depth AS
    SELECT g.representative AS gene_id, d.library_id, SUM(d.aligned_bases / c.length) AS depth
    FROM gene g JOIN contig c USING (contig_id) JOIN contig_dna_counts d USING (contig_id)
    GROUP BY g.representative, d.library_id;

-- Expression over gene copy number, for an RNA and a DNA library of the same sample
CREATE VIEW expression_ratio AS
    SELECT e.gene_id, r.sample_id, e.library_id AS rna_library, d.library_id AS dna_library,
           e.tpm, d.depth AS dna_depth, e.tpm / NULLIF(d.depth, 0) AS ratio
    FROM gene_expression e JOIN library r ON r.library_id = e.library_id
    JOIN library l ON l.sample_id = r.sample_id AND l.molecule = 'dna'
    JOIN gene_dna_depth d ON d.gene_id = e.gene_id AND d.library_id = l.library_id;

-- How much of each DNA library the final MAGs hold
CREATE VIEW library_capture AS
    SELECT l.library_id, l.sample_id, l.primary_mapped, COALESCE(SUM(f.reads), 0) AS reads_on_final_mags,
           COALESCE(SUM(f.reads), 0) / NULLIF(l.primary_mapped, 0) AS fraction_on_final_mags
    FROM library l LEFT JOIN final_mag_dna f USING (library_id)
    WHERE l.molecule = 'dna'
    GROUP BY l.library_id, l.sample_id, l.primary_mapped;
"""

CHECKS = [
    ("mag.representative", "SELECT count(*) FROM mag m WHERE m.representative IS NOT NULL "
                           "AND m.representative NOT IN (SELECT mag_id FROM mag)"),
    ("mag.fine_representative", "SELECT count(*) FROM mag m WHERE m.fine_representative IS NOT NULL "
                                "AND m.fine_representative NOT IN (SELECT mag_id FROM mag)"),
    ("gene.representative", "SELECT count(*) FROM gene g WHERE g.representative NOT IN (SELECT gene_id FROM gene)"),
    ("catalog_sequence.gene_id",
     "SELECT count(*) FROM (SELECT DISTINCT representative FROM gene) r "
     "WHERE r.representative NOT IN (SELECT gene_id FROM catalog_sequence)"),
]


def parse_args():
    parser = argparse.ArgumentParser(description="load every result of the run into one DuckDB database.")
    parser.add_argument("--assemblies", required=True, help="assemblies.tsv: assembly, kind, samples, source.")
    parser.add_argument("--contig-names", nargs="+", required=True, help="<name>.contig_names.tsv of every assembly.")
    parser.add_argument("--contig2mag", nargs="+", required=True, help="<name>.contig2mag.tsv of every assembly.")
    parser.add_argument("--mags-summary", required=True, help="mags_summary.tsv.")
    parser.add_argument("--gff", nargs="+", required=True, help="predicted genes of every assembly, <name>.gff.")
    parser.add_argument("--gene-clusters", required=True, help="gene_clusters.tsv.")
    parser.add_argument("--catalog-faa", required=True, help="gene_catalog.faa.")
    parser.add_argument("--catalog-fna", required=True, help="gene_catalog.fna.")
    parser.add_argument("--counts", help="directory of the tables written by aggregate_counts.py and link_contigs.py.")
    return parser.parse_args()


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, header, rows):
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join("" if v is None else str(v) for v in row) + "\n")


def fasta(path):
    """yield (id, sequence)"""
    name, seq = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def stage_tables(args, tmp):
    """write the tables that need more than a column mapping as TSVs; return their paths"""
    mag_of = {}
    for path in args.contig2mag:
        for row in read_tsv(path):
            mag_of[row["contig"]] = row["mag_id"]
    contig_rows = []
    for path in sorted(args.contig_names):
        assembly = os.path.basename(path).removesuffix(".contig_names.tsv")
        for row in read_tsv(path):
            contig_rows.append((row["contig"], assembly, row["original_contig"], row["length"], mag_of.get(row["contig"])))
    write_tsv(f"{tmp}/contig.tsv", ["contig_id", "assembly_id", "original_contig", "length", "mag_id"], contig_rows)

    clusters = {row["member"]: (row["representative"], row["step"]) for row in read_tsv(args.gene_clusters)}
    gene_rows = []
    for path in sorted(args.gff):
        with open(path) as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                attributes = dict(a.split("=", 1) for a in fields[8].strip(";").split(";") if "=" in a)
                gene = attributes["ID"]
                if gene not in clusters:
                    sys.exit(f"error: gene '{gene}' is in no cluster of gene_clusters.tsv")
                representative, step = clusters[gene]
                gene_rows.append((gene, fields[0], fields[3], fields[4], fields[6], attributes.get("partial"),
                                  representative, step))
    write_tsv(f"{tmp}/gene.tsv", ["gene_id", "contig_id", "start", "stop", "strand", "partial", "representative", "step"],
              gene_rows)

    nucleotide = dict(fasta(args.catalog_fna))
    write_tsv(f"{tmp}/catalog_sequence.tsv", ["gene_id", "protein", "nucleotide"],
              ((gene, protein, nucleotide[gene]) for gene, protein in fasta(args.catalog_faa)))

    assemblies = read_tsv(args.assemblies)
    links = [(row["assembly"], s) for row in assemblies for s in row["samples"].split(";") if s]
    library_samples = []
    if args.counts and os.path.exists(f"{args.counts}/libraries.tsv"):
        library_samples = [row["sample"] for row in read_tsv(f"{args.counts}/libraries.tsv")]
    samples = sorted({s for _, s in links} | set(library_samples))
    write_tsv(f"{tmp}/sample.tsv", ["sample_id"], ((s,) for s in samples))
    write_tsv(f"{tmp}/assembly_sample.tsv", ["assembly_id", "sample_id"], sorted(set(links)))


def load(con, table, path, select):
    """insert the rows of a TSV, every column read as text and cast by `select`"""
    con.execute(f"INSERT INTO {table} SELECT {select} FROM read_csv(?, delim='\\t', header=true, all_varchar=true, "
                f"quote='', escape='')", [path])


def main():
    args = parse_args()
    tmp = "staged"
    os.makedirs(tmp, exist_ok=True)
    stage_tables(args, tmp)

    con = duckdb.connect("clustermg.duckdb")
    con.execute(SCHEMA)
    load(con, "assembly", args.assemblies, "assembly, NULLIF(kind, ''), source")
    load(con, "sample", f"{tmp}/sample.tsv", "sample_id")
    load(con, "assembly_sample", f"{tmp}/assembly_sample.tsv", "assembly_id, sample_id")
    counts = args.counts
    if counts and os.path.exists(f"{counts}/libraries.tsv"):
        load(con, "library", f"{counts}/libraries.tsv",
             "library, sample, molecule, platform, layout, NULLIF(strandedness, ''), input_reads::BIGINT, "
             "NULLIF(rrna_reads, '')::BIGINT, primary_reads::BIGINT, primary_mapped::BIGINT, "
             "NULLIF(assigned, '')::DOUBLE, NULLIF(unassigned_no_features, '')::DOUBLE")
    load(con, "mag", args.mags_summary,
         "mag_id, assembly, original_bin, n_contigs::INTEGER, genome_size::BIGINT, n50::BIGINT, "
         "completeness::DOUBLE, contamination::DOUBLE, status, NULLIF(secondary_cluster, ''), "
         "NULLIF(representative, ''), NULLIF(drep_score, '')::DOUBLE, NULLIF(fine_cluster, ''), "
         "NULLIF(fine_representative, '')")
    load(con, "contig", f"{tmp}/contig.tsv",
         "contig_id, assembly_id, original_contig, length::BIGINT, NULLIF(mag_id, '')")
    load(con, "gene", f"{tmp}/gene.tsv",
         "gene_id, contig_id, start::BIGINT, stop::BIGINT, strand, NULLIF(partial, ''), representative, step")
    load(con, "catalog_sequence", f"{tmp}/catalog_sequence.tsv", "gene_id, protein, nucleotide")
    if counts and os.path.exists(f"{counts}/libraries.tsv"):
        # BED-style links, 0-based and end-exclusive, stored 1-based and inclusive like the genes
        load(con, "contig_link", f"{counts}/contig_links.tsv",
             "link_id, contig, start::BIGINT + 1, \"end\"::BIGINT, mag_id, identity::DOUBLE, aligned_length::BIGINT")
        load(con, "contig_dna_counts", f"{counts}/contig_dna_counts.tsv",
             "contig, library, reads::BIGINT, aligned_bases::BIGINT, covered_bases::BIGINT")
        load(con, "link_dna_counts", f"{counts}/link_dna_counts.tsv",
             "link_id, library, reads::BIGINT, aligned_bases::BIGINT, covered_bases::BIGINT")
        load(con, "mag_dna_profile", f"{counts}/mag_dna_profile.tsv",
             "mag_id, library, reads::BIGINT, mean::DOUBLE, trimmed_mean::DOUBLE, covered_fraction::DOUBLE")
        load(con, "gene_rna_counts", f"{counts}/gene_rna_counts.tsv", "gene, library, reads::DOUBLE")
    con.execute(VIEWS)

    tables = [row[0] for row in con.execute(
        "SELECT table_name FROM duckdb_tables() ORDER BY table_name").fetchall()]
    broken = [(name, con.execute(query).fetchone()[0]) for name, query in CHECKS]
    with open("database_checks.tsv", "w") as out:
        out.write("check\tvalue\n")
        for table in tables:
            out.write(f"rows:{table}\t{con.execute(f'SELECT count(*) FROM {table}').fetchone()[0]}\n")
        for name, value in broken:
            out.write(f"broken:{name}\t{value}\n")
    if any(value for _, value in broken):
        sys.exit(f"error: links that point nowhere: {', '.join(f'{n} ({v})' for n, v in broken if v)}")

    con.execute("EXPORT DATABASE 'parquet' (FORMAT PARQUET)")
    con.close()


if __name__ == "__main__":
    main()
