# dsamoht/clustermg: Usage

## Introduction

clustermg takes metagenome assemblies and the MAGs binned from them, and returns:

- **one dereplicated set of MAGs**, the best MAG of each genome across all the assemblies;
- **one non-redundant gene catalog** (faa, fna, gff) of every gene of every assembly, in which a
  gene of a final MAG is only ever represented by a gene of a final MAG.

It reads the outputs of [mag-ont](https://github.com/dsamoht/mag-ont) as they are, and the bins of
any other binner.

## Samplesheet input

Pass a comma-separated file with four columns to `--input`, one row per assembly.

```bash
--input '[path to samplesheet file]'
```

```csv title="samplesheet.csv"
name,assembly,bins,checkm2
G1,mag-ont/group_G1/assembly/flye/G1.assembly.fasta,mag-ont/group_G1/binning/final_bins,mag-ont/group_G1/binning/checkm/quality_report.tsv
G2,mag-ont/group_G2/assembly/flye/G2.assembly.fasta,mag-ont/group_G2/binning/final_bins,mag-ont/group_G2/binning/checkm/quality_report.tsv
S3,assemblies/S3.contigs.fa.gz,,
```

| Column     | Description                                                                                                                                                                       |
| ---------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `name`     | Assembly name, unique in the samplesheet. Letters, digits, `.`, `_` and `-` only: it prefixes every contig, gene and MAG id of the assembly.                                     |
| `assembly` | The assembly the bins were made from, FASTA (`.fa`, `.fasta`, `.fna`), gzipped or not. Every contig contributes genes to the catalog, binned or not.                            |
| `bins`     | Optional directory of bin FASTA files (`.fa`, `.fasta`, `.fna`, gzipped or not), one file per MAG. Other files in it are ignored.                                                |
| `checkm2`  | CheckM2 `quality_report.tsv` covering the bins: its `Name` column is the bin file name without extension. Required with `bins`. Leave both empty for an assembly without bins. |

A relative path is looked up in the launch directory first, then in the directory of the samplesheet,
so a samplesheet kept next to its data can be used from anywhere.

Checked before any work is done:

- every bin contig is a contig of the assembly, with the same sequence (case aside);
- no contig is in two bins, and no contig name appears twice in an assembly;
- the assembly is DNA;
- every bin has a row in the CheckM2 report;
- an assembly is listed once.

### Names

Contig names collide across assemblies (`k141_1`, `contig_1`), so every contig is renamed
`<name>_<contig>`, where `<contig>` is the first word of its FASTA header, and every gene is named
after its contig, `<name>_<contig>_<n>`. A MAG is named `<name>_<bin file name>`, unless its file name
already starts with `<name>_`, as mag-ont's `<group>_bin_001.fa` does. `assemblies/<name>.contig_names.tsv`
maps the new contig names back to the original ones.

## What the pipeline does

1. **Gene prediction.** Pyrodigal in metagenomic mode on every assembly.
2. **MAG dereplication.**
   - MAGs below `--mag_min_completeness`, above `--mag_max_contamination` or shorter than
     `--mag_min_length` are left out of the final set (`low_quality`).
   - dRep clusters the others at `--mag_ani` with skani, and keeps the best-scoring MAG of each
     cluster (`final`) over the others (`redundant`).
   - With a single MAG past the filter there is nothing to compare, and it is final on its own.
3. **Gene catalog.** Every gene of every assembly is dereplicated with MMseqs2, in three steps that
   make the genes of the final MAGs win:
   1. the final-MAG genes are clustered among themselves;
   2. every other gene, from unbinned contigs and from `redundant` or `low_quality` MAGs, is searched
      against the representatives of step 1, and joins its best hit;
   3. the genes no final-MAG representative took are clustered among themselves.

   The catalog holds the representative of every cluster. The same thresholds, `--gene_min_seq_id`
   and `--gene_coverage`, apply to every step.

### The forcing rule

In every cluster that contains a final-MAG gene, the representative is a final-MAG gene. A gene of a
final MAG is therefore never represented by an unbinned gene or a gene of a discarded MAG. The
pipeline checks this before publishing the catalog and stops if it does not hold.

Two consequences to keep in mind:

- **Final-MAG genes collapse into each other.** Two final MAGs that share a near-identical gene (a
  mobile element, an antibiotic resistance gene) get one catalog entry for it, from one of the two.
  The gene of the other MAG is listed in `gene_clusters.tsv`, under its representative, but is not in
  the catalog faa/fna/gff. `gene_catalog_stats.tsv` reports how many final-MAG contigs still have at
  least one gene of their own in the catalog.
- **A truncated MAG gene can win over a complete one.** A gene cut by the end of a contig can be the
  representative of its cluster. A full-length copy elsewhere that it covers for less than
  `--gene_coverage` of its length cannot join it, and forms a cluster of its own: a little redundancy,
  more frequent with fragmented short-read MAGs than with long-read ones.

The representatives cannot be swapped after a plain clustering of every gene instead: MMseqs2 only
guarantees that each member meets the thresholds against its own representative, so another member
put in its place could leave the rest of the cluster below them.

## Running the pipeline

```bash
nextflow run dsamoht/clustermg -profile docker --input ./samplesheet.csv --outdir ./results
```

The work directory and Nextflow cache are written to the current directory.

### Updating the pipeline

```bash
nextflow pull dsamoht/clustermg
```

### Reproducibility

Specify the pipeline version with `-r`, and record it in your methods:

```bash
nextflow run dsamoht/clustermg -r 2.0.0 -profile docker --input ./samplesheet.csv --outdir ./results
```

## Core Nextflow arguments

### `-profile`

Use one of the container profiles rather than installing the tools yourself. Profiles are
comma-separated and order matters, the last one wins.

- `docker`, `singularity`, `apptainer`, `podman`, `shifter`, `charliecloud`, `wave`
- `conda` / `mamba` — only if none of the container engines are available
- `test` — a minimal dataset that exercises the whole pipeline

### Resource limits

Memory and time requests grow with each retry, and are capped so that no task asks for more than a
machine has: 256 GB, 32 CPUs and 72 h by default. Set `process.resourceLimits` in a config passed with
`-c` to match your largest node.

### `-resume`

Restart a previous run, reusing the tasks whose inputs have not changed.

## Pipeline specific options

| Option                    | Default | Description                                                                                                                      |
| ------------------------- | ------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `--mag_ani`               | `0.95`  | ANI above which two MAGs are the same genome (dRep `-sa`). 0.95 keeps one MAG per species, 0.99 keeps strains apart.             |
| `--mag_min_completeness`  | `50`    | Minimum CheckM2 completeness (%) of a final MAG.                                                                                 |
| `--mag_max_contamination` | `10`    | Maximum CheckM2 contamination (%) of a final MAG.                                                                                |
| `--mag_min_length`        | `50000` | Minimum genome size (bp) of a final MAG (dRep `-l`).                                                                             |
| `--gene_min_seq_id`       | `0.95`  | Minimum amino acid identity between a gene and its representative (MMseqs2 `--min-seq-id`).                                      |
| `--gene_coverage`         | `0.9`   | Minimum alignment coverage (MMseqs2 `-c`), of the sequence `--gene_cov_mode` names.                                              |
| `--gene_cov_mode`         | `1`     | `1`: the member must be covered by its representative, so fragments join full-length genes. `0`: both sequences. `2`: the representative. |

The MAG thresholds are applied by clustermg before dRep and given to dRep as well, so that the
statuses of `mags_summary.tsv` and dRep agree.
