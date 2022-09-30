################################################################################
#
# Investigating biases in microbiome sequencing data: 
#   a meta-analysis of the mock community MSA-2002
#
################################################################################

# ------------------------------------------------------------------------------
# Set up the workspace
# ------------------------------------------------------------------------------

R.version.string # R version 4.0.2

# Package installation
#install.packages("devtools")
#devtools::install_github("benjjneb/dada2", ref="v1.16") 
#install.packages("tidyverse")
#install.packages("reshape2")    
#install.packages("vegan")

# Load required packages 
library(dada2); packageVersion("dada2")
library(tidyverse)
library(reshape2)
library(vegan)

# Set path to raw data
path <- "C:/Users/rauerlui/PhD/Mock meta-analysis/DADA2"
path_all <- list.dirs(path, recursive = F)
length(path_all)

# ------------------------------------------------------------------------------
# DADA2 processing of raw sequencing data
# ------------------------------------------------------------------------------

# Define file name patterns per study
pattern_F_list <- list(rep("1.fastq", 3), "R1_001.fastq", "1.fastq") 
pattern_R_list <- list(rep("2.fastq", 3), "R2_001.fastq", "2.fastq")

# Define truncation, trimming, and paired parameters per study
truncLen_list <- list (c(205, 157), c(163, 111), c(245), c(288, 250), 
                       c(290, 250)) 
trimLeft_list <- list(c(0, 0), c(0, 0), c(0), c(17, 20), c(17, 21)) 
paired_list <- list(TRUE, TRUE, FALSE, TRUE, TRUE)

# Create output lists 
sample.names_list <- list()
seqtab_list <- list()
seqtab.nochim_list <- list()
track_list <- list()

# Loop for joint DADA2 processing of all 17 samples of 5 studies
for (i in 1:length(path_all)) {
  # Set input directory 
  path_i <- path_all[i]
  print(paste0("i = ", i))
  
  # Set DADA2 parameters
  pattern_F_i <- pattern_F_list[[i]]
  pattern_R_i <- pattern_R_list[[i]]
  truncLen_i <- truncLen_list[[i]]
  trimLeft_i <- trimLeft_list[[i]]
  paired_i <- paired_list[[i]]
  
  # Load foward and reverse files 
  fnFs <- sort(list.files(path_i, pattern = pattern_F_i, full.names = TRUE,
                          recursive = TRUE))
  if(paired_i) {
    fnRs <- sort(list.files(path_i, pattern = pattern_R_i, full.names = TRUE,
                            recursive = TRUE))
  }
  
  # Set sample names
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  print(paste0("sample.names = ", sample.names)) 
  
  # Create and save quality profiles
  pdf(paste0(path, "/Quality_plot_F_", i, ".pdf"))
  print(plotQualityProfile(fnFs) + 
          geom_vline(aes(xintercept = trimLeft_i[1]), colour = "blue", 
                     alpha = 0.5) +
          geom_vline(aes(xintercept = truncLen_i[1]), colour = "blue", 
                     alpha = 0.5))
  dev.off()
  if(paired_i) {
    pdf(paste0(path, "/Quality_plot_R_", i, ".pdf"))
    print(plotQualityProfile(fnRs) +
            geom_vline(aes(xintercept = trimLeft_i[2]), colour = "blue", 
                       alpha = 0.5) +
            geom_vline(aes(xintercept = truncLen_i[2]), colour = "blue", 
                       alpha = 0.5))
    dev.off()
  }
  
  # Define names and directory of filtered files in a subdirectory
  filtFs <- file.path(path_i, "filtered", 
                      paste0(sample.names, "_F_filt.fastq.gz"))
  names(filtFs) <- sample.names
  if(paired_i) {
    filtRs <- file.path(path_i, "filtered", 
                        paste0(sample.names, "_R_filt.fastq.gz"))
    names(filtRs) <- sample.names
  }
  
  # Quality filtering
  if (i == 2) {
    out <- filterAndTrim(
      fwd = fnFs, filt = filtFs, 
      if(paired_i) rev = fnRs else rev = NULL, 
      if(paired_i) filt.rev = filtRs else filt.rev = NULL, 
      if(paired_i) maxEE = c(2,2) else maxEE = 2,
      if(paired_i) truncQ = c(2,2) else truncQ = 2,
      truncLen = truncLen_i, trimLeft = trimLeft_i, maxN = 0, rm.phix = TRUE,
      compress = TRUE, multithread = FALSE, verbose = TRUE,
      # To fix study 2 (with different number of sequences in forw./rev. files)
      matchIDs = TRUE) 
  } else {
    out <- filterAndTrim(
      fwd = fnFs, filt = filtFs, 
      if(paired_i) rev = fnRs else rev = NULL, 
      if(paired_i) filt.rev = filtRs else filt.rev = NULL, 
      if(paired_i) maxEE = c(2,2) else maxEE = 2,
      if(paired_i) truncQ = c(2,2) else truncQ = 2,
      truncLen = truncLen_i, trimLeft = trimLeft_i, maxN = 0, rm.phix = TRUE,
      compress = TRUE, multithread = FALSE, verbose = TRUE)
  }
  
  # Learn error rates 
  errF <- learnErrors(filtFs, multithread = FALSE, verbose = 1)
  if(paired_i) {
    errR <- learnErrors(filtRs, multithread = FALSE, verbose = 1)
  }
  
  # Create and save error plots
  pdf(paste0("/Error_plot_F_", i, ".pdf"))
  print(plotErrors(errF, nominalQ = TRUE))
  dev.off()
  if(paired_i) {
    pdf(paste0("/Error_plot_R_", i, ".pdf"))
    print(plotErrors(errR, nominalQ = TRUE))
    dev.off()
  }
  
  # Sample inference 
  dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
  if(paired_i) {
    dadaRs <- dada(filtRs, err = errR, multithread = FALSE)
  }

  # Merge paired reads
  if(paired_i) {
    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
  }
  
  # Create sequence table
  if(paired_i) {
    seqtab <- makeSequenceTable(mergers) 
  } else {
    seqtab <- makeSequenceTable(dadaFs)
  }

  # Remove chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", 
                                      multithread = TRUE, verbose = TRUE)
  
  # Track the loss of reads through the pipeline
  getN <- function(x) sum(getUniques(x))
  # Determine the number of samples for following if condition 
  out_len <- nrow(out)
  # Build track for varying number of samples 
  if (out_len == 1 & paired_i) { # 1 sample, paired data
   track <- cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers), 
                  sum(seqtab.nochim))
  } else if (out_len == 1) { # 1 sample, not paired data
    track <- cbind(out, getN(dadaFs), NA, sum(seqtab), sum(seqtab.nochim))
  } else { # several samples, paired data
    track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
                   sapply(mergers, getN), rowSums(seqtab.nochim))
  }
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                       "nonchim")
  rownames(track) <- sample.names
  
  # Define rownames of output if there is only one sample
  if(out_len == 1) {
    row.names(seqtab) <- sample.names
    row.names(seqtab.nochim) <- sample.names
  }
    
  # Save output
  seqtab_list[[i]] <- seqtab
  seqtab.nochim_list[[i]] <- seqtab.nochim
  track_list[[i]] <- track
  sample.names_list[[i]] <- sample.names
  
  # Prevent R from filtering samples of different studies together
  rm(fnFs, fnRs, filtFs, filtRs)
}

# Combine results for all studies
res_all <- Reduce(function(x, y) 
  merge(x = x, y = y, by = "Sequence",  all = TRUE),
  lapply(seqtab.nochim_list, function(x) 
    rownames_to_column(as.data.frame(t(x)), "Sequence")))
# Change NA to zero in merged results
res_all[is.na(res_all)] <- 0

# Save output
###save(list = ls(), file = "R_objects_DADA2_2022-09-02.RData")
###load(paste0(path, "/R_objects_DADA2_2022-09-02.RData"), verbose = TRUE)

# Assign taxonomic annotation

#get row of sequences of res_all only
res_all_seq <- res_all$Sequence
#View(res_all_seq)
#typeof(res_all_seq)
#as.vector(res_all_seq)
#unlist(res_all_seq, recursive = TRUE, use.names = TRUE)
#as.vector(res_all_seq)

Sys.time() # 4min
set.seed(1)
taxa <- assignTaxonomy(
  res_all$Sequence, paste0(path, "/silva_nr_v132_train_set.fa.gz"), 
  multithread = FALSE)
Sys.time()

# Check for changes in taxon names
table(grepl(pattern = "Phocaeic", taxa)) # not present
table(grepl(pattern = "Bacteroides", taxa)) # 23 ASVs
table(grepl(pattern = "Schaalia", taxa)) # not present
table(grepl(pattern = "Actinomyces", taxa)) # 5 ASVs

# ------------------------------------------------------------------------------
# Taxonomy analysis
# ------------------------------------------------------------------------------

# Rename samples to original names
names(res_all) <- c("Sequence", "ERR6293201", "ERR5463625", "ERR3446055", 
                    "PC4D0", "PC4D1", "PC4D2", "PC4D3", "PC4D4", 
                    "SRR11797017", "SRR11797018", "SRR11797019", "SRR11797020", 
                    "SRR11797021", "SRR11797022", "SRR11797023", "SRR11797025", 
                    "SRR11797026")

# Merge taxa with combined results
taxa_res_all <- merge(res_all, taxa, by.x = "Sequence", by.y = 0, all = TRUE)

# Summarise counts per genus 
res_all_genus <- taxa_res_all %>%
  #Gruppieren
  group_by(Genus) %>%
  #Summe pro Gruppe berechnen
  summarise_if(is.numeric, sum)

# Define expected genera
exp_genera <- c("Acinetobacter", "Schaalia", "Actinomyces", "Bacillus",
                "Bacteroides", "Phocaeicola", "Bifidobacterium", "Clostridium",
                "Cutibacterium", "Deinococcus", "Enterococcus", "Escherichia",
                "Helicobacter", "Lactobacillus", "Neisseria", "Porphyromonas",
                "Pseudomonas", "Rhodobacter", "Staphylococcus", "Streptococcus")
res_all_genus["Genus2"] <- res_all_genus$Genus
# Define all non-expected genera as "Other"
res_all_genus$Genus2[!apply(
  sapply(exp_genera, function(x) grepl(x, res_all_genus$Genus)),
  1, any)] <- "Other"
# List all unique genera
sort(unique(res_all_genus$Genus2)) # 2 types of Clostridium
# Remove detailed information from Clostridium
res_all_genus$Genus2 <- 
  str_replace(res_all_genus$Genus2, "_sensu_stricto_.*", "")
# Remove "Shigella" information from Escherichia
res_all_genus$Genus2 <- 
  str_replace(res_all_genus$Genus2, "\\/Shigella", "")

# Taxonomy plot
taxonomy_plot <- res_all_genus %>%
  # summarise others
  group_by(Genus2) %>%
  summarise_if(is.numeric, sum) %>% #View
  # Add expected column
  add_column(Expected = c(rep(0.05, 13), 0, rep(0.05, 3), rep(0.1, 2))) %>%
  # Create long format
  melt(.) %>% 
  # Add empty column for spacing
  #add_row(Genus2 = "Other", variable = "", value = NA) %>%
  # Order taxonomy levels and samples
  mutate(Genus2 = factor(Genus2, levels = c(
    exp_genera[exp_genera %in% .$Genus2], "Other"))) %>% 
  #mutate(variable = factor(variable, levels = c("Expected", "", colnames(res_all)))) %>% 
  ggplot(., aes(x = variable, y = value, fill = Genus2)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual("Genus", 
                    values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', 
                               '#f58231', '#911eb4', '#46f0f0', '#f032e6', 
                               '#bcf60c', '#fabebe', '#008080', '#e6beff', 
                               '#9a6324', '#fffac8', '#800000', '#aaffc3', 
                               '#808000', '#ffd8b1', #'#000075', '#808080', 
                               #'#ffffff', 
                               '#000000')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0, 0)) + 
  ylab("Relative abundance")
taxonomy_plot

# ------------------------------------------------------------------------------
# Beta diversity analysis
# ------------------------------------------------------------------------------

# Normalise to 100%
res_all_genus_scaled <- res_all_genus[, 1:18]
res_all_genus_scaled[, 2:18] <- sweep(
  res_all_genus_scaled[, 2:18],
  2, colSums(res_all_genus_scaled[, 2:18]), "/")

# Run non-metric multidimensional scaling (nMDS) analysis
set.seed(1)
res_nMDS <- metaMDS(vegdist(t(res_all_genus_scaled[, 2:18])), 
                      try = 100, k = 2)

# Define data for plotting nMDS results
res_nMDS_meta <- data.frame(
  res_nMDS$points,
  Study_ID = c("1", "2", "3", rep("4", 5), rep("5", 9)),
  Extraction_kit = c("A", "B", "A", rep("C", 5), rep("A", 3), rep("D", 3),
                     rep("E", 3)),
  V_region = c("V4", "V3", "V4", rep("V3-V4", 14)),
  Use_beads = c(rep("Yes", 14), rep("No", 3))
)

# Plot nMDS results by Study ID
ggplot(res_nMDS_meta, 
       aes(x = MDS1, y = MDS2, colour = Study_ID)) +
  geom_point() + stat_ellipse()
# Plot nMDS results by Extraction kit
ggplot(res_nMDS_meta, 
       aes(x = MDS1, y = MDS2, colour = Extraction_kit)) +
  geom_point() + stat_ellipse()
# Plot nMDS results by 16S variable region
ggplot(res_nMDS_meta, 
       aes(x = MDS1, y = MDS2, colour = V_region)) +
  geom_point() + stat_ellipse()
# Plot nMDS results by use of beads
ggplot(res_nMDS_meta, 
       aes(x = MDS1, y = MDS2, colour = Use_beads)) +
  geom_point() + stat_ellipse()
