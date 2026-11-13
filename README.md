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

Dereplicate MAGs and genes from metagenomes

     Github: https://github.com/dsamoht/clustermg
     Version: v2.0.0
```

## tl;dr
```bash
nextflow run dsamoht/clustermg --help # nextflow pulls the repo automatically
```

## How it works
Inputs :
  1) assemblies (fasta) 
  2) corresponding MAGs (fasta or "contig2bin" .txt file). 
  3) MAGs statistics from CheckM2  
  4) (optional) DNA reads 
  5) (optional) RNA reads 

clustermg build: 
- one final set of MAGs, dereplicated with dRep after a CheckM2 quality filter;
- one gene catalog of every gene of every assembly, dereplicated with MMseqs2, in which the genes of the final MAGs are only ever represented by a final MAG gene;
- (optional) features (MAGs, genes) abundance/expression level matrices when reads libraries are provided; 
- (experimental) one DuckDB database and a local application to ease results exploration.

## pipeline overview

```mermaid
flowchart LR
    input[/"assemblies, bins,<br/>CheckM2 reports"/] --> PREPARE_ASSEMBLY["PREPARE_ASSEMBLY<br/>rename contigs"]
    PREPARE_ASSEMBLY --> PYRODIGAL["PYRODIGAL<br/>predict genes"]
    reads[/"read libraries<br/>(optional)"/]

    subgraph MAG_DEREPLICATION
        FILTER_MAGS["FILTER_MAGS<br/>CheckM2 quality"] --> DREP["DREP<br/>95% ANI"] --> MAG_LEVELS["MAG_LEVELS<br/>98% ANI groups"] --> SUMMARIZE_MAGS["SUMMARIZE_MAGS<br/>final MAGs"]
    end

    subgraph GENE_CATALOG
        MMSEQS_EASYCLUSTER_MAG["MMSEQS_EASYCLUSTER_MAG<br/>1. final-MAG genes"] --> MMSEQS_EASYSEARCH["MMSEQS_EASYSEARCH<br/>2. other genes"] --> MMSEQS_EASYCLUSTER_OTHER["MMSEQS_EASYCLUSTER_OTHER<br/>3. genes left over"] --> BUILD_CATALOG["BUILD_CATALOG<br/>gene catalog"]
    end

    subgraph READ_MAPPING
        PREPARE_REFERENCE["PREPARE_REFERENCE<br/>all contigs"] --> MINIMAP2_ALIGN_CONTIGS["MINIMAP2_ALIGN_CONTIGS<br/>unbinned vs MAG contigs"] --> LINK_CONTIGS
        SORTMERNA["SORTMERNA<br/>remove rRNA"]
        subgraph map ["map each library once, onto all assemblies"]
            BWAMEM2_MEM["BWAMEM2_MEM<br/>Illumina"]
            MINIMAP2_ALIGN["MINIMAP2_ALIGN<br/>ONT, PacBio"]
        end
        subgraph count ["count"]
            SAMTOOLS_FLAGSTAT["SAMTOOLS_FLAGSTAT<br/>reads"]
            SAMTOOLS_BEDCOV["SAMTOOLS_BEDCOV<br/>DNA on contigs"]
            COVERM_GENOME["COVERM_GENOME<br/>DNA on MAGs"]
            SUBREAD_FEATURECOUNTS["SUBREAD_FEATURECOUNTS<br/>RNA on genes"]
        end
        PREPARE_REFERENCE --> map
        SORTMERNA --> BWAMEM2_MEM
        map --> count
        LINK_CONTIGS --> SAMTOOLS_BEDCOV
        count --> AGGREGATE_COUNTS["AGGREGATE_COUNTS<br/>sum per cluster"]
    end

    PREPARE_ASSEMBLY --> FILTER_MAGS
    PYRODIGAL --> MMSEQS_EASYCLUSTER_MAG
    SUMMARIZE_MAGS --> MMSEQS_EASYCLUSTER_MAG
    PYRODIGAL --> PREPARE_REFERENCE
    SUMMARIZE_MAGS --> PREPARE_REFERENCE
    reads -. DNA .-> map
    reads -. RNA .-> SORTMERNA
    BUILD_CATALOG --> AGGREGATE_COUNTS
    BUILD_CATALOG --> BUILD_DATABASE
    AGGREGATE_COUNTS --> BUILD_DATABASE
    BUILD_DATABASE --> db[("clustermg.duckdb")]
```

## software dependencies

- Nextflow
- a container engine (i.e. docker/singularity) or conda

## documentation

Full documentation: [usage](./docs/usage.md).  
Pipeline outputs: [output](./docs/output.md).
