"""Tests of the helper scripts in bin/, run the way the modules call them.

    pip install pytest
    pytest tests/bin
"""

import csv
import gzip
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[2] / "bin"


def run(script, *args, cwd):
    return subprocess.run(
        [sys.executable, str(BIN / script), *map(str, args)],
        cwd=cwd, capture_output=True, text=True,
    )


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def headers(path):
    return [line[1:].split()[0] for line in Path(path).read_text().splitlines() if line.startswith(">")]


CHECKM_HEADER = ("Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\t"
                 "Coding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\t"
                 "Total_Coding_Sequences\tTotal_Contigs\tMax_Contig_Length\tAdditional_Notes\n")


def checkm_row(name, completeness, contamination):
    return f"{name}\t{completeness}\t{contamination}\tm\t11\t0.9\t4\t300\t8\t0.5\t2\t2\t4\tNone\n"


# ------------------------------------------------------------------------ prepare_assembly

ASSEMBLY = ">c1 flag=1 len=8\nACGTACGT\n>c2 len=4\nGGCC\n>c3\nTTTTAAAA\n"


def prepare(tmp_path, assembly=ASSEMBLY, bins=None, checkm=None, name="S1", assembly_name="asm.fa"):
    write(tmp_path / assembly_name, assembly)
    args = ["--name", name, "--assembly", assembly_name]
    if bins is not None:
        for file_name, text in bins.items():
            write(tmp_path / "bins" / file_name, text)
        (tmp_path / "bins").mkdir(exist_ok=True)
        write(tmp_path / "checkm2.tsv", CHECKM_HEADER + (checkm or ""))
        args += ["--bins", "bins", "--checkm2", "checkm2.tsv"]
    return run("prepare_assembly.py", *args, cwd=tmp_path)


def test_contigs_and_mags_are_renamed(tmp_path):
    result = prepare(tmp_path,
                     bins={"bin.1.fa": ">c1 flag=1 len=8\nACGT\nACGT\n", "S1_bin_002.fasta": ">c3\nTTTTAAAA\n"},
                     checkm=checkm_row("bin.1", 90, 1) + checkm_row("S1_bin_002", 55.5, 4.2))
    assert result.returncode == 0, result.stderr

    assert headers(tmp_path / "S1.contigs.fa") == ["S1_c1", "S1_c2", "S1_c3"]
    assert read_tsv(tmp_path / "S1.contig_names.tsv")[1] == {"contig": "S1_c2", "original_contig": "c2", "length": "4"}
    # a bin already named after its assembly keeps its name, any other gets the prefix
    assert headers(tmp_path / "mags" / "S1_bin.1.fa") == ["S1_c1"]
    assert headers(tmp_path / "mags" / "S1_bin_002.fa") == ["S1_c3"]
    # bins are read in byte order of their file names: "S1_bin_002.fasta" < "bin.1.fa"
    assert [(r["contig"], r["mag_id"]) for r in read_tsv(tmp_path / "S1.contig2mag.tsv")] == [
        ("S1_c3", "S1_bin_002"), ("S1_c1", "S1_bin.1")]
    info = {r["mag_id"]: r for r in read_tsv(tmp_path / "S1.mag_info.tsv")}
    assert info["S1_bin.1"]["original_bin"] == "bin.1.fa"
    assert info["S1_bin.1"]["genome_size"] == "8"
    assert info["S1_bin_002"]["completeness"] == "55.5"
    assert info["S1_bin_002"]["contamination"] == "4.2"


def test_gzipped_assembly_and_bins(tmp_path):
    with gzip.open(tmp_path / "asm.fa.gz", "wt") as f:
        f.write(ASSEMBLY)
    (tmp_path / "bins").mkdir()
    with gzip.open(tmp_path / "bins" / "bin.1.fa.gz", "wt") as f:
        f.write(">c2\nGGCC\n")
    write(tmp_path / "checkm2.tsv", CHECKM_HEADER + checkm_row("bin.1", 90, 1))

    result = run("prepare_assembly.py", "--name", "S1", "--assembly", "asm.fa.gz",
                 "--bins", "bins", "--checkm2", "checkm2.tsv", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert headers(tmp_path / "mags" / "S1_bin.1.fa") == ["S1_c2"]


def test_assembly_without_bins(tmp_path):
    result = prepare(tmp_path)
    assert result.returncode == 0, result.stderr
    assert headers(tmp_path / "S1.contigs.fa") == ["S1_c1", "S1_c2", "S1_c3"]
    assert read_tsv(tmp_path / "S1.contig2mag.tsv") == []
    assert read_tsv(tmp_path / "S1.mag_info.tsv") == []


def test_bin_contig_case_does_not_matter(tmp_path):
    result = prepare(tmp_path, bins={"b.fa": ">c1\nacgtacgt\n"}, checkm=checkm_row("b", 90, 1))
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("bins, checkm, message", [
    ({"b.fa": ">c9\nACGT\n"}, checkm_row("b", 90, 1), "missing from assembly"),
    ({"b.fa": ">c1\nACGTACGA\n"}, checkm_row("b", 90, 1), "differs from the same contig"),
    ({"b.fa": ">c1\nACGTACGT\n", "d.fa": ">c1\nACGTACGT\n"}, checkm_row("b", 90, 1) + checkm_row("d", 90, 1),
     "is in two bins"),
    ({"b.fa": ">c1\nACGTACGT\n>c1\nACGTACGT\n"}, checkm_row("b", 90, 1), "is listed twice in bin"),
    ({"b.fa": ">c1\nACGTACGT\n"}, checkm_row("other", 90, 1), "no CheckM2 row"),
    ({"notes.txt": "not a bin\n"}, "", "no FASTA file"),
    ({"b.fa": ">c1\nACGTACGT\n", "b.fasta": ">c2\nGGCC\n"}, checkm_row("b", 90, 1), "would both become MAG"),
])
def test_bins_that_do_not_match_their_assembly_are_refused(tmp_path, bins, checkm, message):
    result = prepare(tmp_path, bins=bins, checkm=checkm)
    assert result.returncode != 0
    assert message in result.stderr


def test_protein_assembly_is_refused(tmp_path):
    result = prepare(tmp_path, assembly=">p1\nMKVLAAGIVGLLLAQEFP\n")
    assert result.returncode != 0
    assert "is not DNA" in result.stderr


def test_repeated_contig_name_is_refused(tmp_path):
    result = prepare(tmp_path, assembly=">c1\nACGT\n>c1 again\nGGCC\n")
    assert result.returncode != 0
    assert "appears twice" in result.stderr


# ----------------------------------------------------------------------------- filter_mags

MAG_INFO_HEADER = "mag_id\tassembly\toriginal_bin\tn_contigs\tgenome_size\tn50\tcompleteness\tcontamination\n"


def mag_info(*rows):
    return MAG_INFO_HEADER + "".join("\t".join(map(str, r)) + "\n" for r in rows)


def filter_mags(tmp_path, *tables):
    paths = [write(tmp_path / f"t{i}.mag_info.tsv", t).name for i, t in enumerate(tables)]
    return run("filter_mags.py", "--min-completeness", 50, "--max-contamination", 10, "--min-length", 1000,
               *paths, cwd=tmp_path)


def test_filter_keeps_mags_at_the_thresholds(tmp_path):
    result = filter_mags(tmp_path,
                         mag_info(("S1_a", "S1", "a.fa", 1, 1000, 1000, 50, 10),   # exactly at every threshold
                                  ("S1_b", "S1", "b.fa", 1, 999, 999, 90, 1),     # too short
                                  ("S1_c", "S1", "c.fa", 1, 5000, 5000, 49.9, 1)),
                         mag_info(("S2_d", "S2", "d.fa", 1, 5000, 5000, 95, 10.1),
                                  ("S2_e", "S2", "e.fa", 1, 5000, 5000, "", "")))  # not assessed
    assert result.returncode == 0, result.stderr

    assert {r["mag_id"]: r["passes_filter"] for r in read_tsv(tmp_path / "mag_quality.tsv")} == {
        "S1_a": "yes", "S1_b": "no", "S1_c": "no", "S2_d": "no", "S2_e": "no"}
    assert (tmp_path / "passing_mags.txt").read_text() == "S1_a.fa\n"
    assert (tmp_path / "genome_info.csv").read_text() == "genome,completeness,contamination\nS1_a.fa,50,10\n"


def test_filter_refuses_a_mag_id_from_two_assemblies(tmp_path):
    result = filter_mags(tmp_path,
                         mag_info(("X_bin", "X", "bin.fa", 1, 5000, 5000, 90, 1)),
                         mag_info(("X_bin", "X_b", "X_bin.fa", 1, 5000, 5000, 90, 1)))
    assert result.returncode != 0
    assert "comes from assemblies" in result.stderr


# -------------------------------------------------------------------------- summarize_mags

QUALITY_HEADER = MAG_INFO_HEADER.rstrip("\n") + "\tpasses_filter\n"


def summarize(tmp_path, passing, failing=(), drep=None):
    rows = [(m, "S1", f"{m}.fa", 1, 5000, 5000, 90, 1, "yes") for m in passing]
    rows += [(m, "S1", f"{m}.fa", 1, 5000, 5000, 10, 1, "no") for m in failing]
    write(tmp_path / "mag_quality.tsv", QUALITY_HEADER + "".join("\t".join(map(str, r)) + "\n" for r in rows))
    write(tmp_path / "S1.contig2mag.tsv",
          "contig\tmag_id\n" + "".join(f"{m}_ctg1\t{m}\n{m}_ctg2\t{m}\n" for m in [*passing, *failing]))
    mags = [write(tmp_path / "mags" / f"{m}.fa", f">{m}_ctg1\nACGT\n>{m}_ctg2\nACGT\n") for m in [*passing, *failing]]
    args = ["--quality", "mag_quality.tsv", "--contig2mag", "S1.contig2mag.tsv",
            "--mags", *[f"mags/{p.name}" for p in mags]]
    if drep:
        clusters, scores, winners = drep
        write(tmp_path / "Cdb.csv", "genome,secondary_cluster,threshold,cluster_method,comparison_algorithm,primary_cluster\n"
              + "".join(f"{g}.fa,{c},0.05,average,skani,{c.split('_')[0]}\n" for g, c in clusters.items()))
        write(tmp_path / "Sdb.csv", "genome,score\n" + "".join(f"{g}.fa,{s}\n" for g, s in scores.items()))
        write(tmp_path / "Wdb.csv", "genome,cluster,score\n" + "".join(f"{g}.fa,{c},{scores[g]}\n" for c, g in winners.items()))
        # one finer group per dRep cluster, represented by its winner
        write(tmp_path / "mag_levels.tsv", "mag_id\tsecondary_cluster\tfine_cluster\tfine_representative\n"
              + "".join(f"{g}\t{c}\t{c}.1\t{winners[c]}\n" for g, c in clusters.items()))
        args += ["--cdb", "Cdb.csv", "--sdb", "Sdb.csv", "--wdb", "Wdb.csv", "--levels", "mag_levels.tsv"]
    return run("summarize_mags.py", *args, cwd=tmp_path)


def test_summary_with_drep_tables(tmp_path):
    result = summarize(tmp_path, passing=["A", "A2", "C"], failing=["B"],
                       drep=({"A": "1_1", "A2": "1_1", "C": "2_1"},
                             {"A": 80.5, "A2": 60.1, "C": 70.0},
                             {"1_1": "A", "2_1": "C"}))
    assert result.returncode == 0, result.stderr

    summary = {r["mag_id"]: r for r in read_tsv(tmp_path / "mags_summary.tsv")}
    assert {m: r["status"] for m, r in summary.items()} == {
        "A": "final", "A2": "redundant", "B": "low_quality", "C": "final"}
    assert summary["A2"]["representative"] == "A"
    assert summary["A2"]["secondary_cluster"] == "1_1"
    assert summary["A2"]["drep_score"] == "60.1"
    assert (summary["A2"]["fine_cluster"], summary["A2"]["fine_representative"]) == ("1_1.1", "A")
    assert summary["B"]["representative"] == ""
    assert summary["B"]["fine_representative"] == ""
    assert [(r["contig"], r["mag_id"]) for r in read_tsv(tmp_path / "final_contig2mag.tsv")] == [
        ("A_ctg1", "A"), ("A_ctg2", "A"), ("C_ctg1", "C"), ("C_ctg2", "C")]
    assert sorted(p.name for p in (tmp_path / "final_mags").iterdir()) == ["A.fa", "C.fa"]


def test_summary_without_drep_keeps_the_only_passing_mag(tmp_path):
    result = summarize(tmp_path, passing=["A"], failing=["B"])
    assert result.returncode == 0, result.stderr
    summary = {r["mag_id"]: r for r in read_tsv(tmp_path / "mags_summary.tsv")}
    assert {m: r["status"] for m, r in summary.items()} == {"A": "final", "B": "low_quality"}
    assert (summary["A"]["fine_cluster"], summary["A"]["fine_representative"]) == ("", "A")
    assert sorted(p.name for p in (tmp_path / "final_mags").iterdir()) == ["A.fa"]


def test_summary_without_any_mag(tmp_path):
    write(tmp_path / "mag_quality.tsv", QUALITY_HEADER)
    write(tmp_path / "S1.contig2mag.tsv", "contig\tmag_id\n")
    result = run("summarize_mags.py", "--quality", "mag_quality.tsv", "--contig2mag", "S1.contig2mag.tsv",
                 cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert read_tsv(tmp_path / "mags_summary.tsv") == []
    assert read_tsv(tmp_path / "final_contig2mag.tsv") == []


def test_summary_refuses_several_passing_mags_without_drep(tmp_path):
    result = summarize(tmp_path, passing=["A", "C"])
    assert result.returncode != 0
    assert "no dRep table" in result.stderr


def test_summary_refuses_drep_tables_that_disagree_with_the_filter(tmp_path):
    result = summarize(tmp_path, passing=["A", "C"],
                       drep=({"A": "1_1"}, {"A": 80.0}, {"1_1": "A"}))
    assert result.returncode != 0
    assert "dRep clustered" in result.stderr


# ------------------------------------------------------------------------------ mag_levels

NDB_HEADER = "reference,querry,ani,alignment_coverage,primary_cluster\n"


def ndb_rows(primary, ani, coverage=100.0):
    """both directions of every pair in `ani`, plus the self comparisons, as dRep writes them"""
    genomes = sorted({g for pair in ani for g in pair})
    rows = [f"{g}.fa,{g}.fa,1.0,100.0,{primary}\n" for g in genomes]
    for (a, b), value in ani.items():
        rows += [f"{a}.fa,{b}.fa,{value},{coverage},{primary}\n", f"{b}.fa,{a}.fa,{value},{coverage},{primary}\n"]
    return "".join(rows)


def levels(tmp_path, ndb, clusters, scores, fine_ani=0.98, threshold="0.050000000000000044"):
    pytest.importorskip("scipy")
    write(tmp_path / "Ndb.csv", NDB_HEADER + ndb)
    write(tmp_path / "Cdb.csv", "genome,secondary_cluster,threshold,cluster_method,comparison_algorithm,primary_cluster\n"
          + "".join(f"{g}.fa,{c},{threshold},average,skani,{c.split('_')[0]}\n" for g, c in clusters.items()))
    write(tmp_path / "Sdb.csv", "genome,score\n" + "".join(f"{g}.fa,{s}\n" for g, s in scores.items()))
    return run("mag_levels.py", "--ndb", "Ndb.csv", "--cdb", "Cdb.csv", "--sdb", "Sdb.csv",
               "--fine-ani", fine_ani, cwd=tmp_path)


def test_finer_groups_nest_in_drep_clusters(tmp_path):
    # A, A2 and A3 are one species (>= 95%); A3 is a strain of its own below 98%. C is alone.
    ani = {("A", "A2"): 0.99, ("A", "A3"): 0.96, ("A2", "A3"): 0.965}
    result = levels(tmp_path, ndb_rows(1, ani) + "C.fa,C.fa,1.0,100.0,2\n",
                    clusters={"A": "1_1", "A2": "1_1", "A3": "1_1", "C": "2_1"},
                    scores={"A": 70, "A2": 80, "A3": 90, "C": 50})
    assert result.returncode == 0, result.stderr
    rows = {r["mag_id"]: r for r in read_tsv(tmp_path / "mag_levels.tsv")}
    # numbered by the score of their representative: A3 (90) alone, then A2 (80) with A
    assert {m: (r["fine_cluster"], r["fine_representative"]) for m, r in rows.items()} == {
        "A": ("1_1.2", "A2"), "A2": ("1_1.2", "A2"), "A3": ("1_1.1", "A3"), "C": ("2_1.1", "C")}


def test_tree_that_does_not_give_back_drep_clusters_is_refused(tmp_path):
    # 96% is one species, but dRep is said to have kept A and A2 apart
    result = levels(tmp_path, ndb_rows(1, {("A", "A2"): 0.96}),
                    clusters={"A": "1_1", "A2": "1_2"}, scores={"A": 70, "A2": 80})
    assert result.returncode != 0
    assert "does not give back dRep's clusters" in result.stderr


def test_low_coverage_comparisons_are_ignored_as_drep_does(tmp_path):
    # 99% ANI over 0.1% of the genome: at or below -nc, the pair is not the same species
    result = levels(tmp_path, ndb_rows(1, {("A", "A2"): 0.99}, coverage=0.1),
                    clusters={"A": "1_1", "A2": "1_2"}, scores={"A": 70, "A2": 80})
    assert result.returncode == 0, result.stderr
    assert {r["mag_id"]: r["fine_cluster"] for r in read_tsv(tmp_path / "mag_levels.tsv")} == {
        "A": "1_1.1", "A2": "1_2.1"}


def test_finer_level_below_drep_threshold_is_refused(tmp_path):
    result = levels(tmp_path, ndb_rows(1, {("A", "A2"): 0.99}), clusters={"A": "1_1", "A2": "1_1"},
                    scores={"A": 70, "A2": 80}, fine_ani=0.9)
    assert result.returncode != 0
    assert "must be at least as strict" in result.stderr


# ------------------------------------------------------------------------- partition_genes

def test_genes_are_split_on_final_mag_contigs(tmp_path):
    write(tmp_path / "S1.faa", ">S1_k141_1_1 # 1 # 9 # 1 # ID=1_1\nMKV\n>S1_k141_1_2 # x\nMAA\nGG\n"
                               ">S1_k141_12_1 # y\nMQQ\n")
    write(tmp_path / "final_contig2mag.tsv", "contig\tmag_id\nS1_k141_1\tS1_bin_001\n")
    result = run("partition_genes.py", "--faa", "S1.faa", "--final-contig2mag", "final_contig2mag.tsv",
                 "--prefix", "S1", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "S1.mag_genes.faa").read_text() == ">S1_k141_1_1 # 1 # 9 # 1 # ID=1_1\nMKV\n>S1_k141_1_2 # x\nMAA\nGG\n"
    # k141_12 is not k141_1, even though its name starts with it
    assert headers(tmp_path / "S1.other_genes.faa") == ["S1_k141_12_1"]


# ---------------------------------------------------------------------------- assign_genes

def hit(query, target, fident, bits):
    return f"{query}\t{target}\t{fident}\t100\t1.0\t1.0\t1e-50\t{bits}\n"


def test_genes_join_their_best_representative(tmp_path):
    write(tmp_path / "a.other_genes.faa", ">g1\nMKV\n>g2 desc\nMAA\n>g3\nMQQ\nMQQ\n")
    write(tmp_path / "b.other_genes.faa", ">g4\nMRR\n")
    write(tmp_path / "hits.m8",
          hit("g1", "r_9", 0.97, 200) + hit("g1", "r_2", 0.99, 250)      # higher bitscore wins
          + hit("g2", "r_5", 0.96, 300) + hit("g2", "r_3", 0.98, 300)    # same bitscore: identity
          + hit("g4", "r_12", 0.99, 300) + hit("g4", "r_1", 0.99, 300))  # full tie: smallest id
    result = run("assign_genes.py", "--hits", "hits.m8", "b.other_genes.faa", "a.other_genes.faa", cwd=tmp_path)
    assert result.returncode == 0, result.stderr

    assert (tmp_path / "assigned.tsv").read_text() == "r_2\tg1\nr_3\tg2\nr_1\tg4\n"
    assert (tmp_path / "unassigned.faa").read_text() == ">g3\nMQQ\nMQQ\n"


def test_without_hits_every_gene_is_unassigned(tmp_path):
    write(tmp_path / "a.other_genes.faa", ">g1\nMKV\n")
    result = run("assign_genes.py", "a.other_genes.faa", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "assigned.tsv").read_text() == ""
    assert (tmp_path / "unassigned.faa").read_text() == ">g1\nMKV\n"


def test_hits_for_unknown_genes_are_refused(tmp_path):
    write(tmp_path / "a.other_genes.faa", ">g1\nMKV\n")
    write(tmp_path / "hits.m8", hit("ghost", "r_1", 0.99, 100))
    result = run("assign_genes.py", "--hits", "hits.m8", "a.other_genes.faa", cwd=tmp_path)
    assert result.returncode != 0
    assert "absent from the inputs" in result.stderr


# --------------------------------------------------------------------------- build_catalog

# S1_c1 is in final MAG F, S1_c3 in the redundant MAG R, S2_c1 in final MAG G; the rest is unbinned
GENES = {
    "S1": {"S1_c1": ["S1_c1_1", "S1_c1_2"], "S1_c2": ["S1_c2_1"], "S1_c3": ["S1_c3_1"]},
    "S2": {"S2_c1": ["S2_c1_1"], "S2_c2": ["S2_c2_1", "S2_c2_2"]},
}
CONTIG2MAG = {"S1": {"S1_c1": "F", "S1_c3": "R"}, "S2": {"S2_c1": "G"}}
STATUS = {"F": "final", "G": "final", "R": "redundant"}

MAG_CLUSTERS = [("S1_c1_1", "S1_c1_1"), ("S1_c1_1", "S2_c1_1"), ("S1_c1_2", "S1_c1_2")]  # G's gene collapses into F's
ASSIGNED = [("S1_c1_1", "S1_c3_1"), ("S1_c1_2", "S2_c2_1")]
OTHER_CLUSTERS = [("S1_c2_1", "S1_c2_1"), ("S1_c2_1", "S2_c2_2")]


def build(tmp_path, mag_clusters=MAG_CLUSTERS, assigned=ASSIGNED, other_clusters=OTHER_CLUSTERS, genes=GENES):
    for assembly, contigs in genes.items():
        faa = fna = gff = ""
        for contig, contig_genes in contigs.items():
            for i, gene in enumerate(contig_genes):
                faa += f">{gene} # {i * 300 + 1} # {i * 300 + 300} # 1 # ID=x\nMKV\n"
                fna += f">{gene} # {i * 300 + 1} # {i * 300 + 300} # 1 # ID=x\nATGAAAGTG\n"
                gff += f"{contig}\tpyrodigal_v3.7.0\tCDS\t{i * 300 + 1}\t{i * 300 + 300}\t9.9\t+\t0\tID={gene};partial=00;\n"
        write(tmp_path / f"{assembly}.faa", faa)
        write(tmp_path / f"{assembly}.fna", fna)
        write(tmp_path / f"{assembly}.gff", "##gff-version  3\n# Sequence Data: seqnum=1\n" + gff)
        write(tmp_path / f"{assembly}.contig2mag.tsv",
              "contig\tmag_id\n" + "".join(f"{c}\t{m}\n" for c, m in CONTIG2MAG.get(assembly, {}).items()))
    write(tmp_path / "mags_summary.tsv", "mag_id\tstatus\n" + "".join(f"{m}\t{s}\n" for m, s in STATUS.items()))

    args = []
    for flag, pairs, name in [("--mag-clusters", mag_clusters, "mag_cluster.tsv"),
                              ("--assigned", assigned, "assigned.tsv"),
                              ("--other-clusters", other_clusters, "other_cluster.tsv")]:
        if pairs is not None:
            write(tmp_path / name, "".join(f"{r}\t{m}\n" for r, m in pairs))
            args += [flag, name]
    args += ["--faa", *[f"{a}.faa" for a in genes], "--fna", *[f"{a}.fna" for a in genes],
             "--gff", *[f"{a}.gff" for a in genes], "--contig2mag", *[f"{a}.contig2mag.tsv" for a in genes],
             "--mags-summary", "mags_summary.tsv"]
    return run("build_catalog.py", *args, cwd=tmp_path)


def test_catalog_holds_one_representative_per_cluster(tmp_path):
    result = build(tmp_path)
    assert result.returncode == 0, result.stderr

    reps = ["S1_c1_1", "S1_c1_2", "S1_c2_1"]
    assert headers(tmp_path / "gene_catalog.faa") == reps
    assert headers(tmp_path / "gene_catalog.fna") == reps
    gff = (tmp_path / "gene_catalog.gff").read_text().splitlines()
    assert gff[0] == "##gff-version 3"
    assert [line.split("\t")[8].split(";")[0] for line in gff[1:]] == [f"ID={r}" for r in reps]

    clusters = {r["member"]: r for r in read_tsv(tmp_path / "gene_clusters.tsv")}
    assert len(clusters) == 7
    assert clusters["S2_c1_1"] == {"representative": "S1_c1_1", "member": "S2_c1_1", "member_assembly": "S2",
                                   "member_contig": "S2_c1", "member_mag": "G", "member_mag_status": "final",
                                   "step": "mag_cluster"}
    assert clusters["S1_c3_1"]["member_mag_status"] == "redundant"
    assert clusters["S1_c3_1"]["step"] == "mag_assigned"
    assert clusters["S2_c2_2"]["member_mag_status"] == "unbinned"
    assert clusters["S2_c2_2"]["step"] == "de_novo"

    stats = {r["metric"]: int(r["value"]) for r in read_tsv(tmp_path / "gene_catalog_stats.tsv")}
    assert stats == {
        "genes": 7,
        "final_mag_genes": 3,
        "representatives": 3,
        "final_mag_representatives": 2,
        "genes_absorbed_by_final_mag_representatives": 2,
        "de_novo_clusters": 1,
        "final_mag_contigs": 2,
        # G's only gene collapsed into F's, so S2_c1 has no gene of its own in the catalog
        "final_mag_contigs_in_catalog_gff": 1,
    }


def test_catalog_without_any_final_mag(tmp_path):
    genes = {"S1": {"S1_c2": ["S1_c2_1", "S1_c2_2"]}}
    result = build(tmp_path, mag_clusters=None, assigned=[], genes=genes,
                   other_clusters=[("S1_c2_1", "S1_c2_1"), ("S1_c2_1", "S1_c2_2")])
    assert result.returncode == 0, result.stderr
    assert headers(tmp_path / "gene_catalog.faa") == ["S1_c2_1"]


def test_final_mag_gene_represented_from_outside_the_final_mags_is_refused(tmp_path):
    # G's gene ends up under an unbinned representative: exactly what the forcing rule forbids
    result = build(tmp_path,
                   mag_clusters=[("S1_c1_1", "S1_c1_1"), ("S1_c1_2", "S1_c1_2")],
                   other_clusters=OTHER_CLUSTERS + [("S1_c2_1", "S2_c1_1")])
    assert result.returncode != 0
    assert "represented by a gene outside the final MAGs" in result.stderr
    assert not (tmp_path / "gene_catalog.faa").exists()


@pytest.mark.parametrize("change, message", [
    (dict(assigned=ASSIGNED + [("S1_c2_1", "S2_c2_2")]), "member of two clusters"),
    (dict(other_clusters=[("S1_c2_1", "S1_c2_1")]), "belong to no cluster"),
    (dict(other_clusters=[("S1_c2_1", "S2_c2_2"), ("S2_c2_2", "S1_c2_1")]), "is itself a member of another cluster"),
    (dict(assigned=ASSIGNED + [("S1_c1_1", "ghost_1")]), "which no assembly predicts"),
])
def test_inconsistent_clusters_are_refused(tmp_path, change, message):
    result = build(tmp_path, **change)
    assert result.returncode != 0
    assert message in result.stderr


def test_gene_predicted_in_two_assemblies_is_refused(tmp_path):
    genes = {"S1": {"X_c1": ["X_c1_1"]}, "S2": {"X_c1": ["X_c1_1"]}}
    result = build(tmp_path, mag_clusters=None, assigned=[], other_clusters=[("X_c1_1", "X_c1_1")], genes=genes)
    assert result.returncode != 0
    assert "is predicted in both" in result.stderr


# ----------------------------------------------------------------------- prepare_reference

def reference_inputs(tmp_path):
    """two assemblies: S1 with a final MAG (c1) and an unbinned contig (c2), S2 with a low-quality MAG (c1)"""
    write(tmp_path / "S1.contigs.fa", ">S1_c1\nACGTACGT\n>S1_c2\nGGCC\n")
    write(tmp_path / "S2.contigs.fa", ">S2_c1\nTTTTAAAA\nCC\n")
    write(tmp_path / "S1.contig2mag.tsv", "contig\tmag_id\nS1_c1\tS1_m\n")
    write(tmp_path / "S2.contig2mag.tsv", "contig\tmag_id\nS2_c1\tS2_low\n")
    write(tmp_path / "mags_summary.tsv", "mag_id\tassembly\tgenome_size\tstatus\trepresentative\tfine_representative\n"
          "S1_m\tS1\t8\tfinal\tS1_m\tS1_m\nS2_low\tS2\t10\tlow_quality\t\t\n")
    write(tmp_path / "S1.gff", "##gff-version 3\n"
          "S1_c1\tpyrodigal\tCDS\t1\t6\t.\t+\t0\tID=S1_c1_1;partial=00\n"
          "S1_c2\tpyrodigal\tCDS\t1\t3\t.\t-\t0\tID=S1_c2_1;partial=10\n")
    write(tmp_path / "S2.gff", "##gff-version 3\nS2_c1\tpyrodigal\tCDS\t2\t9\t.\t+\t0\tID=S2_c1_1;partial=00\n")


def test_reference_holds_every_contig_and_splits_the_link_alignment(tmp_path):
    reference_inputs(tmp_path)
    result = run("prepare_reference.py", "--contigs", "S1.contigs.fa", "S2.contigs.fa", "--gff", "S1.gff", "S2.gff",
                 "--contig2mag", "S1.contig2mag.tsv", "S2.contig2mag.tsv", "--mags-summary", "mags_summary.tsv",
                 cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert headers(tmp_path / "reference.fa") == ["S1_c1", "S1_c2", "S2_c1"]
    # the contigs of passing MAGs are the target; unbinned and low-quality contigs the query
    assert headers(tmp_path / "link_target.fa") == ["S1_c1"]
    assert headers(tmp_path / "link_query.fa") == ["S1_c2", "S2_c1"]
    assert (tmp_path / "contigs.bed").read_text() == "S1_c1\t0\t8\tcontig\nS1_c2\t0\t4\tcontig\nS2_c1\t0\t10\tcontig\n"
    assert read_tsv(tmp_path / "reference_genes.saf")[1] == {
        "GeneID": "S1_c2_1", "Chr": "S1_c2", "Start": "1", "End": "3", "Strand": "-"}


# ---------------------------------------------------------------------------- link_contigs

def paf(query, qlen, qstart, qend, target, matches, block):
    return f"{query}\t{qlen}\t{qstart}\t{qend}\t+\t{target}\t5000\t0\t{qend - qstart}\t{matches}\t{block}\t60\n"


def test_stretches_link_to_their_best_mag(tmp_path):
    write(tmp_path / "contigs.bed", "U1\t0\t3000\tcontig\nU2\t0\t900\tcontig\n")
    write(tmp_path / "S1.contig2mag.tsv", "contig\tmag_id\nA_c1\tA\nB_c1\tB\n")
    write(tmp_path / "links.paf",
          paf("U1", 3000, 0, 1000, "A_c1", 990, 1000)       # 99% to A
          + paf("U1", 3000, 200, 1200, "B_c1", 970, 1000)   # overlaps the better one: dropped
          + paf("U1", 3000, 2000, 2800, "B_c1", 790, 800)   # 98.75% to B, elsewhere on U1
          + paf("U1", 3000, 1500, 1900, "A_c1", 400, 400)   # too short (400 bp)
          + paf("U2", 900, 0, 900, "A_c1", 800, 900))       # 88.9%: below the identity
    result = run("link_contigs.py", "--paf", "links.paf", "--contigs-bed", "contigs.bed",
                 "--contig2mag", "S1.contig2mag.tsv", "--min-identity", 0.95, "--min-length", 500, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    links = read_tsv(tmp_path / "contig_links.tsv")
    assert [(r["link_id"], r["mag_id"], r["identity"]) for r in links] == [
        ("U1:1-1000", "A", "0.9900"), ("U1:2001-2800", "B", "0.9875")]
    assert (tmp_path / "regions.bed").read_text().splitlines()[2:] == [
        "U1\t0\t1000\tU1:1-1000", "U1\t2000\t2800\tU1:2001-2800"]


def test_without_a_passing_mag_the_regions_are_the_contigs(tmp_path):
    write(tmp_path / "contigs.bed", "U1\t0\t3000\tcontig\n")
    write(tmp_path / "S1.contig2mag.tsv", "contig\tmag_id\n")
    result = run("link_contigs.py", "--contigs-bed", "contigs.bed", "--contig2mag", "S1.contig2mag.tsv",
                 "--min-identity", 0.95, "--min-length", 500, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert read_tsv(tmp_path / "contig_links.tsv") == []
    assert (tmp_path / "regions.bed").read_text() == "U1\t0\t3000\tcontig\n"


# ------------------------------------------------------------------------ aggregate_counts

def aggregate_inputs(tmp_path, rna_assigned=4):
    """final MAG A with redundant A2 (another assembly) and a stretch of unbinned U1 linked to A2"""
    write(tmp_path / "libraries.tsv", "library\tsample\tmolecule\tplatform\tlayout\tstrandedness\n"
          "d1\ts1\tdna\tont\tsingle\t\nr1\ts1\trna\tillumina\tpaired\treverse\n")
    for library in ("d1", "r1"):
        write(tmp_path / f"{library}.flagstat.tsv", "100\t0\tprimary\n90\t0\tprimary mapped\n")
    write(tmp_path / "d1.bedcov.tsv",
          "A_c1\t0\t100\tcontig\t500\t80\t5\n"       # final MAG A
          "A2_c1\t0\t200\tcontig\t300\t60\t3\n"      # redundant A2, represented by A
          "U1\t0\t1000\tcontig\t900\t500\t9\n"       # unbinned contig...
          "U1\t0\t400\tU1:1-400\t200\t100\t2\n"      # ...whose first 400 bp copy A2
          "E\t0\t50\tcontig\t0\t0\t0\n")             # nothing: no row
    write(tmp_path / "contig_links.tsv", "link_id\tcontig\tstart\tend\tmag_id\tidentity\taligned_length\n"
          "U1:1-400\tU1\t0\t400\tA2\t0.9800\t400\n")
    write(tmp_path / "S.contig2mag.tsv", "contig\tmag_id\nA_c1\tA\nA2_c1\tA2\n")
    write(tmp_path / "mags_summary.tsv",
          "mag_id\tassembly\tgenome_size\tstatus\trepresentative\tfine_representative\n"
          "A\tS1\t100\tfinal\tA\tA\nA2\tS2\t200\tredundant\tA\tA2\n")
    write(tmp_path / "gene_clusters.tsv", "representative\tmember\nA_c1_1\tA_c1_1\nA_c1_1\tU1_1\nU1_2\tU1_2\n")
    write(tmp_path / "reference_genes.saf", "GeneID\tChr\tStart\tEnd\tStrand\n"
          "A_c1_1\tA_c1\t1\t100\t+\nU1_1\tU1\t1\t300\t+\nU1_2\tU1\t401\t900\t-\n")
    # a fragment split three ways is written 0.33 each: 0.01 short, within featureCounts' rounding
    write(tmp_path / "r1.featureCounts.tsv", "# Program:featureCounts v2.1.1\n"
          "Geneid\tChr\tStart\tEnd\tStrand\tLength\tr1.bam\n"
          "A_c1_1\tA_c1\t1\t100\t+\t100\t1.33\nU1_1\tU1\t1\t300\t+\t300\t1.33\nU1_2\tU1\t401\t900\t-\t500\t1.33\n")
    write(tmp_path / "r1.featureCounts.tsv.summary",
          f"Status\tr1.bam\nAssigned\t{rna_assigned}\nUnassigned_NoFeatures\t6\n")
    write(tmp_path / "coverm.tsv", "Genome\td1 Read Count\td1 Mean\td1 Trimmed Mean\td1 Covered Fraction\n"
          "A\t5\t5.0\t4.5\t0.8\nA2\t3\t1.5\t1.2\t0.3\n")
    # SortMeRNA was given 130 reads, 30 of them rRNA: the 100 of the flagstat reached the mapper
    write(tmp_path / "r1.sortmerna.log", " Results:\n    Total reads = 130\n"
          "    Total reads passing E-value threshold = 30 (23.08)\n")
    return run("aggregate_counts.py", "--libraries", "libraries.tsv", "--flagstat", "d1.flagstat.tsv", "r1.flagstat.tsv",
               "--bedcov", "d1.bedcov.tsv", "--featurecounts", "r1.featureCounts.tsv", "--coverm", "coverm.tsv",
               "--sortmerna-log", "r1.sortmerna.log",
               "--contig-links", "contig_links.tsv", "--contig2mag", "S.contig2mag.tsv",
               "--mags-summary", "mags_summary.tsv", "--gene-clusters", "gene_clusters.tsv",
               "--saf", "reference_genes.saf", cwd=tmp_path)


def test_counts_add_up_through_the_links(tmp_path):
    result = aggregate_inputs(tmp_path)
    assert result.returncode == 0, result.stderr
    assert len(read_tsv(tmp_path / "contig_dna_counts.tsv")) == 3
    assert read_tsv(tmp_path / "link_dna_counts.tsv") == [
        {"link_id": "U1:1-400", "library": "d1", "reads": "2", "aligned_bases": "200", "covered_bases": "100"}]
    # final MAG A: its contigs, A2's, and the stretch of U1 linked to A2; depth over A's 100 bp
    assert read_tsv(tmp_path / "final_mag_dna.tsv") == [
        {"mag_id": "A", "library": "d1", "reads": "10", "aligned_bases": "1000", "mean_depth": "10"}]
    # at the finer level A2 is a group of its own, with its linked stretch
    assert {r["mag_id"]: (r["reads"], r["mean_depth"]) for r in read_tsv(tmp_path / "fine_group_dna.tsv")} == {
        "A": ("5", "5"), "A2": ("5", "2.5")}
    assert read_tsv(tmp_path / "mag_dna_profile.tsv")[0]["trimmed_mean"] == "4.5"
    # RNA: U1_1 is a copy of A_c1_1, so the cluster has 3 fragments
    expression = {r["representative"]: r for r in read_tsv(tmp_path / "gene_expression.tsv")}
    assert expression["A_c1_1"]["reads"] == "2.66"
    assert round(sum(float(r["tpm"]) for r in expression.values())) == 1_000_000
    libraries = {r["library"]: r for r in read_tsv(tmp_path / "libraries.tsv")}
    assert (libraries["r1"]["assigned"], libraries["r1"]["unassigned_no_features"]) == ("4", "6")
    assert libraries["d1"]["primary_mapped"] == "90"
    assert (libraries["r1"]["input_reads"], libraries["r1"]["rrna_reads"]) == ("130", "30")
    assert (libraries["d1"]["input_reads"], libraries["d1"]["rrna_reads"]) == ("100", "")


def test_gene_counts_that_disagree_with_featurecounts_are_refused(tmp_path):
    result = aggregate_inputs(tmp_path, rna_assigned=5)
    assert result.returncode != 0
    assert "featureCounts assigned 5 fragments" in result.stderr


# -------------------------------------------------------------------------- build_database

def test_database_holds_every_table_and_derives_the_rest(tmp_path):
    pytest.importorskip("duckdb")
    aggregate_inputs(tmp_path)
    counts = tmp_path / "counts"
    counts.mkdir()
    for name in ("libraries", "contig_links", "contig_dna_counts", "link_dna_counts", "mag_dna_profile",
                 "gene_rna_counts"):
        (counts / f"{name}.tsv").write_text((tmp_path / f"{name}.tsv").read_text())
    write(tmp_path / "assemblies.tsv", "assembly\tkind\tsamples\tsource\nS1\tsingle\ts1\ta.fa\nS2\tco\ts1;s2\tb.fa\n")
    write(tmp_path / "S1.contig_names.tsv", "contig\toriginal_contig\tlength\nA_c1\tc1\t100\nU1\tu1\t1000\nE\te\t50\n")
    write(tmp_path / "S2.contig_names.tsv", "contig\toriginal_contig\tlength\nA2_c1\tc1\t200\n")
    write(tmp_path / "mags_full.tsv",
          "mag_id\tassembly\toriginal_bin\tn_contigs\tgenome_size\tn50\tcompleteness\tcontamination\tstatus\t"
          "secondary_cluster\trepresentative\tdrep_score\tfine_cluster\tfine_representative\n"
          "A\tS1\ta.fa\t1\t100\t100\t90\t1\tfinal\t1_1\tA\t80\t1_1.1\tA\n"
          "A2\tS2\ta2.fa\t1\t200\t200\t80\t2\tredundant\t1_1\tA\t60\t1_1.2\tA2\n")
    write(tmp_path / "S1.gff", "##gff-version 3\n"
          "A_c1\tp\tCDS\t1\t100\t.\t+\t0\tID=A_c1_1;partial=00\n"
          "U1\tp\tCDS\t1\t300\t.\t+\t0\tID=U1_1;partial=00\nU1\tp\tCDS\t401\t900\t.\t-\t0\tID=U1_2;partial=01\n")
    write(tmp_path / "S2.gff", "##gff-version 3\n")
    write(tmp_path / "clusters.tsv", "representative\tmember\tstep\nA_c1_1\tA_c1_1\tmag_cluster\n"
          "A_c1_1\tU1_1\tmag_assigned\nU1_2\tU1_2\tde_novo\n")
    write(tmp_path / "catalog.faa", ">A_c1_1\nMK*\n>U1_2\nMV*\n")
    write(tmp_path / "catalog.fna", ">A_c1_1\nATGAAA\n>U1_2\nATGGTT\n")
    result = run("build_database.py", "--assemblies", "assemblies.tsv",
                 "--contig-names", "S1.contig_names.tsv", "S2.contig_names.tsv",
                 "--contig2mag", "S.contig2mag.tsv", "--mags-summary", "mags_full.tsv", "--gff", "S1.gff", "S2.gff",
                 "--gene-clusters", "clusters.tsv", "--catalog-faa", "catalog.faa", "--catalog-fna", "catalog.fna",
                 "--counts", "counts", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    checks = {r["check"]: r["value"] for r in read_tsv(tmp_path / "database_checks.tsv")}
    assert all(v == "0" for k, v in checks.items() if k.startswith("broken:"))
    assert checks["rows:gene"] == "3" and checks["rows:sample"] == "2" and checks["rows:assembly_sample"] == "3"

    import duckdb
    con = duckdb.connect(str(tmp_path / "clustermg.duckdb"), read_only=True)
    # the same sums as aggregate_counts.py, from the raw tables
    assert con.execute("SELECT mag_id, reads, mean_depth FROM final_mag_dna").fetchall() == [("A", 10, 10.0)]
    assert {g: round(r, 2) for g, r in con.execute("SELECT gene_id, reads FROM gene_expression").fetchall()} == {
        "A_c1_1": 2.66, "U1_2": 1.33}
    # U1_1 sits in the stretch of U1 linked to A2, so its fragments count for final MAG A
    assert [(m, round(r, 2)) for m, r in con.execute("SELECT mag_id, reads FROM mag_expression").fetchall()] == [("A", 2.66)]
    assert con.execute("SELECT start, stop FROM contig_link").fetchall() == [(1, 400)]
    assert (tmp_path / "parquet" / "schema.sql").exists()
