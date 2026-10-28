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
    assert read_tsv(tmp_path / "S1.contig_names.tsv")[1] == {"contig": "S1_c2", "original_contig": "c2"}
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
        args += ["--cdb", "Cdb.csv", "--sdb", "Sdb.csv", "--wdb", "Wdb.csv"]
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
    assert summary["B"]["representative"] == ""
    assert [(r["contig"], r["mag_id"]) for r in read_tsv(tmp_path / "final_contig2mag.tsv")] == [
        ("A_ctg1", "A"), ("A_ctg2", "A"), ("C_ctg1", "C"), ("C_ctg2", "C")]
    assert sorted(p.name for p in (tmp_path / "final_mags").iterdir()) == ["A.fa", "C.fa"]


def test_summary_without_drep_keeps_the_only_passing_mag(tmp_path):
    result = summarize(tmp_path, passing=["A"], failing=["B"])
    assert result.returncode == 0, result.stderr
    assert {r["mag_id"]: r["status"] for r in read_tsv(tmp_path / "mags_summary.tsv")} == {
        "A": "final", "B": "low_quality"}
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
