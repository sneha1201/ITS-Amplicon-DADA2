# ITS-Amplicon-DADA2

# ITS DADA2 Pipeline

This repository contains an end-to-end ITS amplicon analysis pipeline using **DADA2** in R.

---

## **Overview**

The pipeline performs the following steps:

1. **Cutadapt Primer Trimming** – trims ITS primers from raw FASTQ reads.
2. **Quality Control** – generates quality profiles for each sample.
3. **Filtering and Trimming** – removes low-quality reads.
4. **Error Learning & Denoising** – models and corrects sequencing errors.
5. **Dereplication & Sequence Inference** – identifies unique sequences (ASVs).
6. **Merging Paired Reads** – merges forward and reverse reads.
7. **Chimera Removal** – removes chimeric sequences.
8. **Taxonomy Assignment** – assigns taxonomy using an ITS reference database (e.g., UNITE).

---

## **Requirements**

- R ≥ 4.0
- Packages: `dada2`, `phyloseq`, `DECIPHER`, `ggplot2`, `Biostrings`
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) installed and in your PATH
- ITS reference database (FASTA) for taxonomy assignment

---

## **Usage**

```bash
Rscript ITS_dada2_pipeline.R <path_to_fastq> <path_to_taxonomy_files> <output_rdata_filename> <output_folder> <cutadapt_path>

Arguments:

path_to_fastq – folder containing raw FASTQ files.

path_to_taxonomy_files – folder containing ITS reference FASTA files.

output_rdata_filename – name for the saved RData object.

output_folder – folder to store results (quality plots, tables, etc.).

cutadapt_path – full path to the cutadapt executable.

Outputs Formats

Quality profiles: quality_profiles/

Filtered FASTQ: filtered_fastq/

Error plots: error_plots/

Sequence tables & ASV counts: seq_tables/

Taxonomy assignments: taxonomy_assignments/

Read tracking table: Number_of_reads_at_each_step.csv

Notes

Raw FASTQ files must end with _R1_001.fastq.gz and _R2_001.fastq.gz.

The script supports multithreading for faster processing.

Use an appropriate ITS reference database compatible with your samples (e.g., UNITE).
