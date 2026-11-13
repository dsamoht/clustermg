# dsamoht/clustermg: Usage

## Introduction

clustermg takes metagenome assemblies and the MAGs binned from them, and returns:

- **one dereplicated set of MAGs**, the best MAG of each genome across all the assemblies;
- **one non-redundant gene catalog** (faa, fna, gff) of every gene of every assembly, in which a
  gene of a final MAG is only ever represented by a gene of a final MAG;
- with read libraries, **the reads of every library counted once**, on every contig and gene, and
  added up into final MAGs and catalog genes (see [Read libraries](#read-libraries));
- **one DuckDB database** holding all of it (see [The database](#the-database)).

Single-sample assemblies and co-assemblies of the same samples can be given together: every genome
and gene comes out once, and every duplicate is kept as a row that points to it.

It reads the outputs of [mag-ont](https://github.com/dsamoht/mag-ont) as they are, and the bins of
any other binner.

## Samplesheet input

Pass a comma-separated file with five columns to `--input`, one row per assembly.

```bash
--input '[path to samplesheet file]'
```
```csv title="samplesheet.csv"
name,assembly,bins,checkm2,samples
co,mag-ont/group_co/assembly/flye/co.assembly.fasta,mag-ont/group_co/binning/final_bins,mag-ont/group_co/binning/checkm/quality_report.tsv,s1;s2
s1,mag-ont/group_s1/assembly/flye/s1.assembly.fasta,mag-ont/group_s1/binning/final_bins,mag-ont/group_s1/binning/checkm/quality_report.tsv,s1
s2,mag-ont/group_s2/assembly/flye/s2.assembly.fasta,mag-ont/group_s2/binning/final_bins,mag-ont/group_s2/binning/checkm/quality_report.tsv,s2
s3,assemblies/s3.contigs.fa.gz,,,s3
```

| Column     | Description                                                                                                                                                                       |
| ---------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `name`     | Assembly name, unique in the samplesheet. Letters, digits, `.`, `_` and `-` only: it prefixes every contig, gene and MAG id of the assembly.                                     |
| `assembly` | The assembly the bins were made from, FASTA (`.fa`, `.fasta`, `.fna`), gzipped or not. Every contig contributes genes to the catalog, binned or not.                            |
| `bins`     | Optional directory of bin FASTA files (`.fa`, `.fasta`, `.fna`, gzipped or not), one file per MAG. Other files in it are ignored.                                                |
| `checkm2`  | CheckM2 `quality_report.tsv` covering the bins: its `Name` column is the bin file name without extension. Required with `bins`. Leave both empty for an assembly without bins. |
| `samples`  | Optional `;`-separated ids of the samples the assembly was built from: one for a single-sample assembly, several for a co-assembly. Provenance only, for the database: it changes no result. |

mag-ont builds both kinds in one run when a sample is listed both in a group of its own and in a
co-assembly group; give each mag-ont group as one row, with its sample ids in `samples`.

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

## Read libraries

Pass a second comma-separated file to `--reads` to count reads, one row per library:

```csv title="reads.csv"
library,sample,molecule,platform,fastq_1,fastq_2,strandedness
s1_dna,s1,dna,ont,reads/s1.ont.fastq.gz,,
s2_dna,s2,dna,ont,reads/s2.ont.fastq.gz,,
s1_rna,s1,rna,illumina,reads/s1_rna_R1.fastq.gz,reads/s1_rna_R2.fastq.gz,reverse
```

| Column         | Description                                                                                                                       |
| -------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| `library`      | Library id, unique. Letters, digits, `.`, `_` and `-` only.                                                                       |
| `sample`       | The sample the library comes from, as in the `samples` of the assemblies. A DNA and an RNA library of one sample share it.      |
| `molecule`     | `dna`, counted on contigs and MAGs, or `rna`, counted on genes.                                                                   |
| `platform`     | `ont`, `pacbio` (taken as HiFi) or `illumina`. RNA must be `illumina`.                                                           |
| `fastq_1`      | FastQ, gzipped or not.                                                                                                           |
| `fastq_2`      | Mates of `fastq_1`, Illumina only.                                                                                               |
| `strandedness` | RNA only, required: `unstranded`, `forward` (read 1 on the transcript's strand) or `reverse` (read 1 antisense, as dUTP kits). |

A FastQ file can be given once only: its reads would otherwise be counted twice. A sample listed in
the `samples` of an assembly without a library gives a warning.

### rRNA removal

The rRNA reads of every RNA library are removed with SortMeRNA before it is mapped: both mates go
when either aligns to an rRNA database. Left in, rRNA reads, often most of a metatranscriptome, would
make mapping slower, and those landing on a predicted gene would be counted as its expression.
`--ribo_database_manifest` lists the databases, one FASTA path or URL per line; by default
SortMeRNA's own (bacterial, archaeal and eukaryotic 16S/18S, 23S/28S, 5S and 5.8S), downloaded when
the run starts. `--skip_rrna_removal` maps the RNA as it is, for reads already filtered. clustermg
does not look for rRNA genes on the contigs: no annotation is run.

### One mapping per library, onto every assembly

Each library is mapped **once** onto the contigs of **all** the assemblies together: with bwa-mem2
for Illumina reads, DNA and RNA alike (prokaryotic genes are not spliced), with minimap2 for long
reads. A read lands on the copy it matches best, the strain's own where several assemblies hold
one, and on one of them at random where the copies are identical. Each read sits on one copy only,
so the copies of a genome or a gene add up without counting it twice:

- **DNA** is counted on every contig: reads, aligned bases and bases at depth 1 or more
  (`samtools bedcov`), from the primary alignment of each read only, so that a chimeric long read
  is counted once. A final MAG holds its whole dRep cluster: the contigs of the final MAG and of
  its `redundant` MAGs, plus the stretches of unbinned or low-quality contigs that align to one of
  them at `--link_min_identity` over `--link_min_length` (a copy of the genome no bin held). Its
  mean depth is the aligned bases of all of them over the final MAG's genome size.
- **RNA** is counted on every gene with featureCounts, strand-aware: one count per fragment, split
  between the genes it overlaps (an operon), so that a library's gene counts add up to its assigned
  fragments. Each gene's count goes to its catalog representative; TPM is computed on the
  representative's length.

Mapping each library separately onto every assembly it helped build would not do: the same read
would be counted once per assembly, by a factor that changes from gene to gene. mag-ont's binning
BAMs are such mappings, which is why the reads are mapped again here.

What one mapping onto every assembly cannot do:

- **Tell identical copies apart.** Where two assemblies hold the same sequence, the read is placed
  on one at random: the counts of each version only mean something where the versions differ. The
  sums are not affected.
- **Keep strains apart in the sums.** A final MAG adds up every strain of its 95% cluster, a finer
  group every strain of its 98% group. The strain detail stays in the counts of each copy.
- **Give a whole cluster a trimmed mean or a covered fraction.** Those are computed per MAG, on its
  own contigs (CoverM); a final MAG's reads and mean depth cover its whole cluster.
- **Count a genome no assembly holds.** Its reads stay unmapped.

The reference holds every assembly, so it is larger than the unique content of the run: about two
to three times with single-sample assemblies and co-assemblies of the same samples. Plan the memory
of the indexes accordingly:

- **bwa-mem2** needs 28 bytes per base to index the reference, about 140 GB for 5 Gb of contigs;
  the task asks for that much, within `process.resourceLimits`.
- **minimap2** indexes a reference larger than its batch size (`-I`, 8 GB) in parts, and merges
  the parts back at mapping.

## The database

`database/clustermg.duckdb` holds every result, each fact once:

- **tables**: `assembly`, `sample`, `assembly_sample`, `library`, `contig`, `mag`, `gene`,
  `catalog_sequence`, `contig_link`, and the raw counts `contig_dna_counts`, `link_dna_counts`,
  `mag_dna_profile` and `gene_rna_counts`;
- **views**, everything that is derived: `final_mag`, `catalog_gene`, `gene_context` (gene →
  contig → assembly → MAG → finer group → final MAG), `final_mag_dna`, `fine_group_dna`,
  `gene_expression`, `mag_expression`, `gene_dna_depth`, `expression_ratio` (RNA TPM over the DNA
  depth of the same sample) and `library_capture` (share of each DNA library on the final MAGs).

clustermg does no annotation. Annotations made outside it (GTDB-Tk on the final MAGs, functions on
the catalog representatives) join on the same ids: annotate the catalog once, and
`gene.representative` carries it to every gene, contig and MAG. `database/parquet/` holds every table as Parquet, with `schema.sql` and `load.sql` to
rebuild the database in DuckDB.

## What the pipeline does

1. **Gene prediction.** Pyrodigal in metagenomic mode on every assembly.
2. **MAG dereplication.**
   - MAGs below `--mag_min_completeness`, above `--mag_max_contamination` or shorter than
     `--mag_min_length` are left out of the final set (`low_quality`).
   - dRep clusters the others at `--mag_ani` with skani, and keeps the best-scoring MAG of each
     cluster (`final`) over the others (`redundant`), whichever kind of assembly it comes from.
   - The same dRep tree, cut at `--mag_ani_fine`, gives finer groups nested in dRep's clusters,
     reported next to them. The tree is rebuilt from dRep's comparisons, and the run stops if it
     does not give back dRep's own clusters.
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
- `test` — a minimal dataset to run the whole pipeline

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
| `--mag_ani_fine`          | `0.98`  | ANI of the finer groups, cut from dRep's tree. 0.98 is where dRep and inStrain consider genomes distinct for read mapping.        |
| `--mag_min_completeness`  | `50`    | Minimum CheckM2 completeness (%) of a final MAG.                                                                                 |
| `--mag_max_contamination` | `10`    | Maximum CheckM2 contamination (%) of a final MAG.                                                                                |
| `--mag_min_length`        | `50000` | Minimum genome size (bp) of a final MAG (dRep `-l`).                                                                             |
| `--gene_min_seq_id`       | `0.95`  | Minimum amino acid identity between a gene and its representative (MMseqs2 `--min-seq-id`).                                      |
| `--gene_coverage`         | `0.9`   | Minimum alignment coverage (MMseqs2 `-c`), of the sequence `--gene_cov_mode` names.                                              |
| `--gene_cov_mode`         | `1`     | `1`: the member must be covered by its representative, so fragments join full-length genes. `0`: both sequences. `2`: the representative. |
| `--reads`                 | —       | Sheet of read libraries, see [Read libraries](#read-libraries).                                                                  |
| `--ribo_database_manifest` | SortMeRNA's databases | rRNA databases the RNA reads are filtered against, one FASTA path or URL per line. See [rRNA removal](#rrna-removal). |
| `--skip_rrna_removal`     | `false` | Map the RNA libraries without removing their rRNA reads.                                                                         |
| `--link_min_identity`     | `0.98`  | Minimum identity of a stretch of an unbinned or low-quality contig to count as part of the MAG it aligns to.                    |
| `--link_min_length`       | `500`   | Minimum length (bp) of such a stretch.                                                                                           |
| `--save_bams`             | `false` | Publish the BAM of every library.                                                                                                |

The MAG thresholds are applied by clustermg before dRep and given to dRep as well, so that the
statuses of `mags_summary.tsv` and dRep agree.
