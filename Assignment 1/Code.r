# Goal: To find taxa for the unclassified sequences with length >300 bp

# Note: All the required files are saved in the working directory for convenience

# ------------------------------
# Step 1: Load Required Libraries
# ------------------------------
# Load essential R packages for data handling, sequence processing, and phylogenetic tree generation.
library(tidyverse)    # For Data manipulation
library(Biostrings)   # For DNA sequence handling
library(DECIPHER)     # For Multiple sequence alignment
library(ape)          # For Phylogenetic tree construction
library(phangorn)     # For Tree analysis

# ------------------------------
# Step 2: Load Kraken2 Data
# ------------------------------
# Read the Kraken2 classification results file
kraken_data <- read_tsv("barcode01.kraken2.assignments.tsv", col_names = FALSE)

# Assign meaningful column names for clarity
colnames(kraken_data) <- c("Status", "Sequence_ID", "Taxonomy_ID", "Length", "KMer_LCA")

# Display first few rows
head(kraken_data)

# ------------------------------
# Step 3: Filter Unclassified Sequences (with length >300 bp)
# ------------------------------
# Extract sequences labeled as "U" and longer than 300 bp
unclassified_seq <- kraken_data %>%
  filter(Status == "U" & Length > 300)

# Check how many unclassified sequences match the condition
nrow(unclassified_seq)
head(unclassified_seq)

# ------------------------------
# Step 4: Extract the respective Sequences from FASTA File
# ------------------------------
# Read the FASTA file containing all sequences
fasta_file <- readDNAStringSet("barcode01_merged.fasta")

# Ensure Sequence_IDs match between Kraken2 output file and FASTA file

# Removing Extra information from FASTA headers
# Remove extra text after the first space in fasta_file
clean_fasta_file <- gsub(" .*", "", names(fasta_file))
head(clean_fasta_file)

# Extract sequences that match unclassified Sequence_IDs
unclassified_fasta <- fasta_file[clean_fasta_file %in% unclassified_seq$Sequence_ID]

# Save extracted unclassified sequences into a new FASTA file
writeXStringSet(unclassified_fasta, "unclassified_sequences.fasta")

# ------------------------------
# Step 5: Perform BLAST Search
# ------------------------------
# Run BLASTn manually by uploading the "unclassified_sequences.fasta" file
# BLASTn Result:
# All sequences show ~95 percent identity with "Enterococcus faecalis"
# Generating a phylogenetic tree to place the sequences in relation to known species, helping confirm their true belonging to E.faecalis or a different species

# ------------------------------
# Step 6: Phylogenetic Analysis to Confirm Taxonomy
# ------------------------------
# Using Neighbor-Joining (NJ) since the sequences are closely related
# Load unclassified sequences again for phylogenetic analysis
unclassified_fasta_saved <- readDNAStringSet("unclassified_sequences.fasta")

# Perform multiple sequence alignment (MSA)
alignment <- AlignSeqs(unclassified_fasta_saved)

# Convert aligned sequences to phylogenetic format
dna_bin <- as.DNAbin(alignment)

# Compute genetic distances
dna_dist <- dist.dna(dna_bin)

# Build a phylogenetic tree using the Neighbor-Joining (NJ) method
nj_tree <- nj(dna_dist)

# Plot the tree
pdf("Unclassified_sequences_nj_tree.pdf", width = 16, height = 7)
plot(nj_tree, main = "Phylogenetic Tree of Unclassified Sequences", cex = 0.7)
dev.off()

# ------------------------------
# Step 7: Compare Unclassified Sequences to Known Reference Sequences 
# (Here, only one Reference Sequence i.e NZ OVVR01000001.1 Enterococcus faecalis isolate Efs173 NODE 1 length 1018502 cov 29.8739, whole genome shotgun sequence is used due to memory constraints)
# ------------------------------
# Load reference sequences
ref_seq <- readDNAStringSet("E.faecalis_Ref_seq.fasta")

# Combine reference sequences with unclassified sequences
all_seq <- c(unclassified_fasta_saved, ref_seq)

# Perform multiple sequence alignment (MSA) on all sequences
alignment_all <- AlignSeqs(all_seq)

# Convert alignment to phylogenetic format
dna_bin_all <- as.DNAbin(alignment_all)

# Compute genetic distances
dna_dist_all <- dist.dna(dna_bin_all)

# Handle missing values in the distance matrix (if any)
if (sum(is.na(dna_dist_all)) > 0) {
  cat("Warning: Distance matrix contains missing values. Removing problematic sequences...\n")
  dna_dist_all_clean <- na.omit(as.matrix(dna_dist_all))
  dna_dist_all <- as.dist(dna_dist_all_clean)
}

# ------------------------------
# Step 8: Construct Parsimony Tree
# ------------------------------
# Using Parsimony for ref. seq. + unclassified sequences bacause of high diversity and unreliable distance calculations
parsimony_tree <- pratchet(as.phyDat(dna_bin_all))

# Save phylogenetic tree as a PDF file
pdf("parsimony_tree.pdf", width = 16, height = 7)
plot(parsimony_tree, main = "Parsimony Phylogenetic Tree", cex = 0.7)
dev.off()

# Interpretation:
# Strong BLAST hits (~95% identity) shows a close relationship to E.faecalis however, its not a ceertainty.
# The unclassified sequences might belong to a related clade as NJ analysis showed high similarity.
# The unclassified sequences are genetically close but aren't identical to E.faecalis.
# This could indicate a possible novel taxa or variation in strain.

# ------------------------------
# Summary of Findings
# ------------------------------
cat("\nSummary of Findings:\n")
cat("- Extracted unclassified sequences (>300 bp) from Kraken2 output.\n")
cat("- Performed BLASTn search to find taxonomic matches.\n")
cat("- Constructed phylogenetic trees to confirm species identity.\n")
cat("- Unclassified sequences show high similarity (~95%) to Enterococcus faecalis.\n")
cat("- Phylogenetic analysis supports classification as E. faecalis.\n")

cat("\nResults saved in 'Unclassified_sequences_nj_tree.pdf' and 'parsimony_tree.pdf' (phylogenetic tree).\n")
