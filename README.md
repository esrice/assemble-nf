# assemble-nf

A minimalist nextflow pipeline for long-read assembly

## Introduction

This pipeline takes long reads for an arbitrary number of samples, and then assembles and evaluates contigs for each sample. For the moment, it accepts hifi reads only, although I'm going to add nanopore functionality in the near future. Currently, it runs the following steps:

* Assemble reads into contigs with hifiasm
* Calculate super basic stats (N50, length) with assembly-stats
* Get BUSCO scores with compleasm
* Get QV scores with inspector

This is meant for cases where you want to assemble a bunch of genomes to contig level without doing the kind of integration of multiple data types or manual curation you'd do for a reference assembly — e.g., for building a pangenome.

## Installation

You just need nextflow and conda or mamba; the pipeline will take care of creating conda envs with all the required software.

## Running

An optional but highly recommended first step is to download the BUSCOs for your lineage so that you only need to do this once total rather than once per assembly:
```
compleasm download aves
library_path=$PWD/mb_downloads
```

(or use the `-L` argument to specify a path other than `./mb_downloads` for the BUSCOs)

Then, make a sample sheet in csv format with fields `sample_id` and `fastqs`. `sample_id` should be a unique identifier for each sample and `fastqs` a semicolon-delimited list of fastq(.gz) files containing long reads to assemble.

Finally, run nextflow:
```
nextflow run esrice/assemble-nf --samplesheet samples.csv --lineage aves --odbLocation $library_path
```

If you're running it on `coolmuc`, I've included a coolmuc profile, so add `-profile coolmuc` to the nextflow command. If you're running it on a different cluster, you'll have to configure nextflow for your cluster, which is outside the scope of this document, but hopefully the included coolmuc profile is a good start.
