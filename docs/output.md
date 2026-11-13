# dsamoht/clustermg: Output

## Introduction

Every path below is relative to
`--outdir`.

## Pipeline overview

- [MAGs](#mags) — the dereplicated MAG set and the fate of every input MAG
- [Gene catalog](#gene-catalog) — the non-redundant genes and the cluster of every gene
- [Gene predictions](#gene-predictions) — every gene of every assembly
- [Assemblies](#assemblies) — the renamed contigs
- [Read counts](#read-counts) — every library on every contig and gene, with `--reads` only
- [Database](#database) — every result in one DuckDB database
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
    - `fine_cluster`, `fine_representative`: the finer group (`--mag_ani_fine`, 98% by default),
      named `<secondary_cluster>.<n>`, and the MAG that represents it. A MAG dereplicated alone is
      its own representative.
  - `final_contig2mag.tsv`: `contig, mag_id` for every contig of the final MAGs.
  - with `--reads`:
    - `final_mag_dna.tsv`: every final MAG × every DNA library: `mag_id, library, reads,
      aligned_bases, mean_depth`. Its whole dRep cluster: the contigs of the final MAG and of its
      `redundant` MAGs, and the stretches of other contigs linked to them; mean depth is over the
      final MAG's genome size.
    - `fine_group_dna.tsv`: the same for every finer group, named after its representative.
    - `mag_dna_profile.tsv`: every input MAG × every DNA library, on its own contigs only (CoverM):
      `mag_id, library, reads, mean, trimmed_mean, covered_fraction`.
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
  - with `--reads`:
    - `gene_expression.tsv`: catalog gene × RNA library: `representative, library, reads, tpm`, the
      fragments of every copy of the gene added up, TPM on the representative's length. Rows with
      no fragment are left out.
    - `gene_rna_counts.tsv`: every gene copy × RNA library: `gene, library, reads`, one count per
      fragment split between the genes it overlaps. Rows with no fragment are left out.
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
  - `<name>.contig_names.tsv`: `contig, original_contig, length` for every contig of the assembly.
  - `<name>.contig2mag.tsv`: `contig, mag_id` for every contig of every MAG of the assembly, final
    or not.

</details>

### Read counts

<details markdown="1">
<summary>Output files</summary>

- `coverage/`, with `--reads` only
  - `libraries.tsv`: `library, sample, molecule, platform, layout, strandedness`, then the reads of
    the library, mates counted separately: `input_reads` (in its FastQ files), `rrna_reads` (removed
    by SortMeRNA before mapping, RNA only), `primary_reads` and `primary_mapped` (samtools
    flagstat); and for RNA, `assigned` and `unassigned_no_features` (featureCounts fragments;
    `unassigned_no_features` is mostly the tRNA, the rRNA left and the intergenic).
  - `contig_dna_counts.tsv`: every contig × DNA library: `contig, library, reads, aligned_bases,
    covered_bases`. Rows with nothing are left out.
  - `contig_links.tsv`: `link_id, contig, start, end, mag_id, identity, aligned_length`: the
    stretches of unbinned or low-quality contigs counted as part of a MAG (`start` 0-based, `end`
    exclusive; `link_id` is `contig:start+1-end`).
  - `link_dna_counts.tsv`: the same counts as `contig_dna_counts.tsv`, on every link.
- `mapping/`, with `--save_bams`
  - `<library>.bam`, `<library>.bam.bai`: the library mapped onto every assembly.

</details>

How the libraries are mapped and counted, and what that cannot tell: see
[Read libraries](./usage.md#read-libraries).

### Database

<details markdown="1">
<summary>Output files</summary>

- `database/`
  - `clustermg.duckdb`: every result, tables for the facts and views for what is derived. See
    [The database](./usage.md#the-database).
  - `parquet/`: every table as Parquet, with `schema.sql` and `load.sql`.
  - `database_checks.tsv`: rows of every table, and the links inside a table that point nowhere
    (`broken:*`, all 0: the run stops otherwise).

</details>

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `clustermg_software_versions.yml`: versions of every tool that ran.
  - `execution_report_*.html`, `execution_timeline_*.html`, `execution_trace_*.txt`,
    `pipeline_dag_*.html`: Nextflow run reports.

</details>
