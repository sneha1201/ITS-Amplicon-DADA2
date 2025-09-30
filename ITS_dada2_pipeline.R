#!/usr/bin/env Rscript

library(dada2)
library(phyloseq)
library(DECIPHER)
library(ggplot2)
library(Biostrings)

# =======================
# Command line arguments
# =======================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript dada2_pipeline.R <path_to_fastq> <path_to_taxonomy_files> <output_rdata_filename> <output_folder>")
}

path <- args[1]
tax_path <- args[2]
rdata_file <- args[3]
outdir <- args[4]
cutadapt <- args[5]

# =======================
# Create output directory
# =======================
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# =======================
# Primer sequences
# =======================
FWD <- "CTTGGTCATTTAGAGGAAGTAA"  
REV <- "GCTGCGTTCTTCATCGATGC" 
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# =======================
# Input FASTQ files
# =======================
fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# =======================
# Cutadapt trimming
# =======================
cut_dir <- file.path(path, "cutadapt")
if (!dir.exists(cut_dir)) dir.create(cut_dir)

cutFs <- file.path(cut_dir, paste0(sample.names, "_R1_cut.fastq.gz"))
cutRs <- file.path(cut_dir, paste0(sample.names, "_R2_cut.fastq.gz"))

for (i in seq_along(fnFs)) {
  cmd <- paste(
    cutadapt,
    "-g", FWD, "-G", REV,          # forward primers
    "-a", REV.RC, "-A", FWD.RC,    # reverse complements
    "-n", "2",                     # search for both primers
    "-o", cutFs[i], "-p", cutRs[i],
    fnFs[i], fnRs[i]
  )
  system(cmd)
}



# =======================
# Quality plots
# =======================
for (i in seq_along(sample.names)) {
  qprofiles <- tryCatch({
    plotQualityProfile(c(cutFs[i], cutRs[i]))
  }, error = function(e) {
    message("Skipping sample ", sample.names[i], ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(qprofiles)) {
    ggsave(filename = file.path(outdir, paste0("quality_profile_", sample.names[i], ".tiff")),
           plot = qprofiles, height = 10, width = 15, units = "in", dpi = 300)
  }
}

# =======================
# Filter and Trim
# =======================
filtered_dir <- file.path(path, "filtered")
if (!dir.exists(filtered_dir)) dir.create(filtered_dir)

filtFs <- file.path(filtered_dir, paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs <- file.path(filtered_dir, paste0(sample.names, "_R2_filt.fastq.gz"))

out <- filterAndTrim(
  cutFs, filtFs, cutRs, filtRs,
  maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 50, rm.phix = TRUE,
  compress = TRUE, multithread = TRUE
)

# Check filtered reads
out

# Check file sizes
file.info(filtFs)$size
file.info(filtRs)$size

# =======================
# Learn Errors
# =======================
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

ggsave(file.path(outdir, "Errors_plot_Forward_Reads.tiff"), plot = plotErrors(errF, nominalQ = TRUE), height = 12, width = 16, dpi = 300)
ggsave(file.path(outdir, "Errors_plot_Reverse_Reads.tiff"), plot = plotErrors(errR, nominalQ = TRUE), height = 12, width = 16, dpi = 300)

# =======================
# Dereplication
# =======================
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# =======================
# DADA2 inference
# =======================
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# =======================
# Merging
# =======================
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# =======================
# Sequence table + Chimera removal
# =======================
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)

# =======================
# Track reads
# =======================
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.csv(track, file = file.path(outdir, "Number_of_reads_at_each_step.csv"))

# =======================
# Taxonomy assignment
# =======================
taxa <- assignTaxonomy(seqtab.nochim, file.path(tax_path, "sh_general_release_dynamic_all_27.10.2022.fasta"), multithread = TRUE)
# taxa <- addSpecies(taxa, file.path(tax_path, "sh_general_release_dynamic_all_19.02.2025.fasta"))

taxa.print <- taxa
rownames(taxa.print) <- NULL
write.csv(as.data.frame(taxa.print), file = file.path(outdir, "taxonomy_assignments.csv"))

# =======================
# Save everything
# =======================
save(list = ls(), file = file.path(outdir, rdata_file))

