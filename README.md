## RNA-seq processing pipeline

This repository contains a pipeline for the analysis of RNA-seq (and derived variants) data. 

__Motivation__

The aim was to create a standard RNA-seq pipeline that included quality controls, alignment and post-processing.
The pipeline has been developed using [Snakemake](http://snakemake.readthedocs.io/en/stable/), a tool to create reproducible and scalable data analyses.

It was used by Federico for his `intergenic transcription` project, and then generalised to become a tool than could be applied to any set of RNA-seq experiments.

__Requirements__

* [Snakemake](http://snakemake.readthedocs.io/en/stable/) with [Python](https://www.python.org/downloads/) >= 3.6 (`os`, `platform`, `pandas`, `numpy`, `string`, `re` and `itertools` modules)
<!--* [_SRA Toolkit_](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software)-->
* [_FastQC_](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [_Trim Galore!_](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) and its dependency [_Cutadapt_](https://cutadapt.readthedocs.io/en/stable/)
* [_Bowtie2_](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [_STAR_](https://github.com/alexdobin/STAR)
* [_StringTie_](https://ccb.jhu.edu/software/stringtie)

## Setup

Please, before running the pipeline modify the `config.yaml` according to your needs.

### Workflow

__Directed acyclic graph (DAG)__

![](img/dag.svg)

Figure 1: Schematic of the pipeline workflow.

---

### Pipeline breakdown

#### Steps description

* Quality control on the raw data
    - Performed with FastQC
* Trimming of the adapter sequence
    - Performed with Trim Galore!
* Removal of ribosomal and transfer RNAs contaminations
    - Performed with Bowtie2
* Alignment to the genome
    - Performed with STAR
* Generation of bigWig for visualisation
    - Performed using the ENCODE binaries
* Transcriptome assembly
    - Performed with Stringtie

#### [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

#### [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

Trim Galore! is a wrapper script to automate quality and adapter trimming as well as quality control, with some added functionality to remove biased methylation positions for RRBS sequence files (for directional, non-directional (or paired-end) sequencing).

#### [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences.

#### [STAR](https://github.com/alexdobin/STAR)

Basic STAR workflow consists of 2 steps:
* Generating genome indexes files.
In this step user supplied the reference genome sequences (FASTA files) and annotations (GTF file), from which STAR generate genome indexes that are utilized in the 2nd (mapping) step. 
* Mapping reads to the genome.
In this step user supplies the genome files generated in the 1st step, as well as the RNA-seq reads (sequences) in the form of FASTA or FASTQ files. STAR maps the reads to the genome,
and writes several output files, such as alignments (SAM/BAM), mapping summary statistics, splice junctions, unmapped reads, signal (wiggle) tracks etc.

#### [StringTie](https://ccb.jhu.edu/software/stringtie)

StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. It uses a novel network flow algorithm as well as an optional de novo assembly step to assemble and quantitate full-length transcripts representing multiple splice variants for each gene locus.

_NB: Presently, the algorithm performs the assembly only on samples that contain some keywords (i.e., WildType or WT, Untreated, siLuc and Senescence)._

---

### To-Do list

  - PCR duplicates removal with Picard (downstream files are currenly produced from STAR outputs);
  - NET/GRO-seq branch (useful for techniques where splicing is not accounted for);
  - Dynamic retrieval of genome and annotation.
