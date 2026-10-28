# dsamoht/clustermg: Output

## Introduction

This document describes the output produced by the pipeline. Every path below is relative to
`--outdir`.

## Pipeline overview

- [MAGs](#mags) — the dereplicated MAG set and the fate of every input MAG
- [Gene catalog](#gene-catalog) — the non-redundant genes and the cluster of every gene
- [Gene predictions](#gene-predictions) — every gene of every assembly
- [Assemblies](#assemblies) — the renamed contigs
- [Pipeline information](#pipeline-information)

### MAGs

<details markdown="1">
<summary>Output files</summary>

- `mags/final/`
  - `<mag_id>.fa`: the final MAGs, one per dRep cluster. Contig headers are the renamed contigs,
    `<name>_<contig>`, the seqids of `genes/gene_catalog.gff`.
- `mags/`
  - `mags_summary.tsv`: every input MAG, one row each, columns in this order:
    - `mag_id`, `assembly`, `original_bin`: the MAG, the assembly it comes from, and its bin file.
    - `n_contigs`, `genome_size`, `n50`: computed from the bin.
    - `completeness`, `contamination`: from the CheckM2 report.
    - `status`: `final`, `redundant` (represented by another MAG of its dRep cluster) or
      `low_quality` (failed the quality filter, never dereplicated).
    - `secondary_cluster`, `representative`, `drep_score`: the dRep cluster, the final MAG that
      represents it, and its dRep score. Empty for `low_quality` MAGs, and every column but
      `representative` is empty when dRep did not run.
  - `final_contig2mag.tsv`: `contig, mag_id` for every contig of the final MAGs.
- `mags/drep/`
  - `Cdb.csv`, `Ndb.csv`, `Sdb.csv`, `Wdb.csv`, `Widb.csv`, `genomeInformation.csv`: dRep's
    clusters, secondary (skani) comparisons, scores, winners, winner details and genome
    statistics. Only when dRep ran, i.e. when at least two MAGs passed the filter.

</details>

### Gene catalog

<details markdown="1">
<summary>Output files</summary>

- `genes/`
  - `gene_catalog.faa`, `gene_catalog.fna`, `gene_catalog.gff`: the representative of every gene
    cluster, with Pyrodigal's headers and attributes. A gene id is `<name>_<contig>_<n>`.
  - `gene_clusters.tsv`: every gene of every assembly, one row each, columns in this order:
    - `representative`, `member`: the catalog gene that represents the gene, and the gene.
    - `member_assembly`, `member_contig`, `member_mag`: where the gene comes from; `member_mag` is
      empty for an unbinned contig.
    - `member_mag_status`: `final`, `redundant`, `low_quality` (the status of its MAG, see
      `mags_summary.tsv`) or `unbinned`.
    - `step`: `mag_cluster` (a final-MAG gene, clustered with the other final-MAG genes),
      `mag_assigned` (another gene, absorbed by a final-MAG representative) or `de_novo` (clustered
      with the genes no final-MAG representative took).
  - `gene_catalog_stats.tsv`: number of genes, final-MAG genes, representatives, final-MAG
    representatives, genes absorbed by final-MAG representatives and de novo clusters; number of
    final-MAG contigs, and of final-MAG contigs with at least one gene in the catalog gff.

</details>

Every gene of a final MAG has a final-MAG gene as representative: the pipeline stops before
publishing the catalog otherwise. Final-MAG genes can represent each other, so a final-MAG contig
can have genes in `gene_clusters.tsv` that are not in the catalog files. See the forcing rule in
[usage](./usage.md#the-forcing-rule).

### Gene predictions

<details markdown="1">
<summary>Output files</summary>

- `genes/predicted/`
  - `<name>.gff`, `<name>.fna`, `<name>.faa`: Pyrodigal (metagenomic mode) predictions on the
    renamed contigs of each assembly, every gene of `gene_clusters.tsv` included.

</details>

### Assemblies

<details markdown="1">
<summary>Output files</summary>

- `assemblies/`
  - `<name>.contig_names.tsv`: `contig, original_contig` for every contig of the assembly.

</details>

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `clustermg_software_versions.yml`: versions of every tool that ran.
  - `execution_report_*.html`, `execution_timeline_*.html`, `execution_trace_*.txt`,
    `pipeline_dag_*.html`: Nextflow run reports.

</details>
