[![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

```
         __              __
  _____ / /__  __ _____ / /_ ___   _____ ____ ___   ____ _
 / ___// // / / // ___// __// _ \ / ___// __ `__ \ / __ `/
/ /__ / // /_/ /(__  )/ /_ /  __// /   / / / / / // /_/ /
\___//_/ \__,_//____/ \__/ \___//_/   /_/ /_/ /_/ \__, /
                                                 /____/

Dereplicated MAGs and a MAG-anchored non-redundant gene catalog
from metagenome assemblies

     Github: https://github.com/dsamoht/clustermg
     Version: v2.0.0
```

## tl;dr

```bash
nextflow run main.nf --help # if working in the repo
```

or

```bash
nextflow run dsamoht/clustermg --help # nextflow pulls the repo automatically
```

From assemblies and the MAGs binned from them (e.g. [mag-ont](https://github.com/dsamoht/mag-ont)
outputs), clustermg builds:

- one set of MAGs, dereplicated with dRep (95% ANI by default) after a CheckM2 quality filter;
- one gene catalog (faa, fna, gff) of every gene of every assembly, dereplicated with MMseqs2 (95%
  identity, 90% coverage by default), in which the genes of the final MAGs always win: a final-MAG
  gene is only ever represented by a final-MAG gene.

## software dependencies

- Nextflow
- a container engine (i.e. docker/singularity) or conda

## documentation

Full documentation: [usage](./docs/usage.md).  
Pipeline outputs: [output](./docs/output.md).
