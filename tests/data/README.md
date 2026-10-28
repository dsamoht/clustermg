# Test data

Small enough to run in a few minutes, built so that every fate a MAG and a gene can have shows up.
The paths of the sample sheets are relative to this directory, where they are resolved.

| File | What it is |
|---|---|
| `assemblies/minigut1.fa` | A MEGAHIT assembly of a mock gut community (269 contigs, 1.7 Mb) |
| `assemblies/minigut2.fa.gz` | A second, fragmented assembly of the same community (1313 contigs, 0.6 Mb), gzipped |
| `assemblies/minigut1_mut.fa` | Synthetic: the contigs of `minigut1_bin_001` with about 1% of their bases substituted |
| `bins/minigut1/minigut1_bin_001.fa` | The longest contigs of `minigut1.fa`, about 600 kb (16 contigs) |
| `bins/minigut1/minigut1_bin_002.fa` | The next longest, about 300 kb (22 contigs) |
| `bins/minigut1_mut/bin.1.fa` | Every contig of `minigut1_mut.fa` |
| `checkm2/*.tsv` | CheckM2 `quality_report.tsv` layout; completeness and contamination are **made up** |

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

- `samplesheet_test.csv`: the three assemblies, `-profile test`.
- `samplesheet_nobins.csv`: assemblies only, no MAG at all.
- `samplesheet_single.csv`: `minigut1` alone, so a single MAG passes the filter and dRep is skipped.

The bins, the mutated assembly and the CheckM2 reports were derived from `minigut1.fa` once, with a
seeded random generator (seed 20261028, 1% substitution rate, bins filled longest contig first);
`minigut2.fa.gz` was compressed with a zero timestamp so that it is byte-identical when rebuilt.
