# Test data

Small enough to run in a few minutes, built so that every fate a MAG and a gene can have shows up.
The paths of the sample sheets are relative to this directory, where they are resolved.

| File | What it is |
|---|---|
| `assemblies/minigut1.fa` | A MEGAHIT assembly of a mock gut community (269 contigs, 1.7 Mb) |
| `assemblies/minigut2.fa.gz` | A second, fragmented assembly of the same community (1313 contigs, 0.6 Mb), gzipped |
| `assemblies/minigut1_mut.fa` | Synthetic: the contigs of `minigut1_bin_001` with about 1% of their bases substituted |
| `assemblies/minigut1_frag.fa.gz` | Synthetic: 116 pieces of 1 to 6 kb of `minigut1_bin_001`'s contigs, 0.5% substituted, never binned |
| `bins/minigut1/minigut1_bin_001.fa` | The longest contigs of `minigut1.fa`, about 600 kb (16 contigs) |
| `bins/minigut1/minigut1_bin_002.fa` | The next longest, about 300 kb (22 contigs) |
| `bins/minigut1_mut/bin.1.fa` | Every contig of `minigut1_mut.fa` |
| `checkm2/*.tsv` | CheckM2 `quality_report.tsv` layout; completeness and contamination are **made up** |
| `reads/s1_dna.fq.gz`, `s2_dna`, `s3_dna` | Simulated ONT-like reads, 150 per library, from `minigut1`, `minigut1_frag` and `minigut1_mut` |
| `reads/s1_rna_R1.fq.gz`, `_R2` | Simulated Illumina RNA, reverse stranded: 2000 read pairs from 200 genes of `minigut1`, then 500 rRNA pairs |
| `rrna/rrna_test.fasta`, `rrna-db-test.txt` | A test rRNA database and its manifest: the first 20 sequences of SortMeRNA's `silva-bac-16s-id90` and the first 10 of `silva-bac-23s-id98` (v4.3.4, derived from SILVA, CC BY 4.0) |

Completeness and contamination given to each bin:

| Bin | Completeness | Contamination | Expected status |
|---|---|---|---|
| `minigut1_bin_001` | 85 | 2 | `final` |
| `minigut1_bin_002` | 40 | 5 | `low_quality` (below 50% completeness) |
| `minigut1_mut/bin.1` (`minigut1_mut_bin.1`) | 70 | 3 | `redundant`, represented by `minigut1_bin_001`: about 99% ANI to it, lower dRep score |

The bin of `minigut1_mut` is named `bin.1.fa`, not after its assembly, so it also exercises the
MAG id prefix rule. Most genes of the mutated copy are expected to be absorbed by the genes of
`minigut1_bin_001` (step `mag_assigned`); the substitutions that create a stop codon leave a few
truncated genes to be clustered on their own.

Sample sheets:

- `samplesheet_test.csv`: the four assemblies, `-profile test`, with the samples each was built
  from: `minigut1` is a co-assembly of s1 and s2, the others single-sample assemblies.
- `reads_test.csv`: the four read libraries, `-profile test`. Each genome copy is reached a different
  way: s1's reads land on `minigut1_bin_001` itself, s3's on its redundant copy `minigut1_mut_bin.1`,
  and s2's on `minigut1_frag`, whose pieces are linked to `minigut1_bin_001`: all three add up into it.
- `reads_duplicate.csv`: one FastQ for two libraries, which must be refused.
- `reads_rna_ont.csv`: an RNA library of ONT reads, which must be refused.
- `samplesheet_nobins.csv`: assemblies only, no MAG at all.
- `samplesheet_single.csv`: `minigut1` alone, so a single MAG passes the filter and dRep is skipped.

The bins, the mutated assembly and the CheckM2 reports were derived from `minigut1.fa` once, with a
seeded random generator (seed 20261028, 1% substitution rate, bins filled longest contig first);
`minigut2.fa.gz` was compressed with a zero timestamp so that it is byte-identical when rebuilt.

`minigut1_frag.fa.gz` was cut from `minigut1_bin_001.fa` with a seeded generator (seed 20261104):
consecutive pieces of 2 to 6 kb, each kept with probability 0.7 if at least 1 kb long, then 0.5% of
its bases substituted. The reads were simulated with seeded generators too, and gzipped with a zero
timestamp:

- DNA (seeds 1 to 3 for s1 to s3): lengths log-normal around 3 kb (500 bp to 8 kb, at most the
  contig), drawn over the contigs of at least 500 bp weighted by length, 2% substitutions, 1%
  insertions, 1% deletions, random strand.
- RNA (seed 4): 200 genes of `minigut1` (Pyrodigal, at least 450 bp with 30 bp on each side), with
  log-normal expression; fragments of 200 to 400 bp within the gene and its sides; read 1 is the
  first 150 bp of the antisense strand of the transcript and read 2 of its sense strand (dUTP
  libraries, `reverse`), 0.5% substitutions. Then 500 pairs drawn the same way from the sequences
  of `rrna/rrna_test.fasta` (seed 5), named `rrna<n>`: SortMeRNA must remove all 1000 reads, and
  the 2000 pairs left must all land on genes.

The gzip headers of the reads hold no file name and a zero timestamp, so that they are
byte-identical when rebuilt.
