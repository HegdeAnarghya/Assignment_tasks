# **Task 1: Taxonomic Classification of Unclassified Sequences**

## **Overview**

This task aims to identify taxa for unclassified sequences from a **Kraken2 classification report**, specifically those with a sequence length **greater than 300 bp**. The approach involves **filtering**, **BLAST searches**, and **phylogenetic analysis** to infer taxonomic relationships.

## **Workflow Summary**

1.  **Extracted Unclassified Sequences:**
    -   Parsed the Kraken2 output (`barcode01.kraken2.assignments.tsv`) to filter sequences labeled **"U"** (Unclassified) with a length **\>300 bp**.
    -   Retrieved corresponding sequences from `barcode01_merged.fasta`.
2.  **BLAST Search for Taxonomic Identification:**
    -   Conducted a **BLASTn search** on `unclassified_sequences.fasta` against the NCBI database.
    -   Results indicated **~95% identity** with _Enterococcus faecalis_.
3.  **Phylogenetic Analysis for Confirmation:**
    -   Used **Multiple Sequence Alignment (MSA)** and **Neighbor-Joining (NJ) method** to construct a phylogenetic tree.
    -   Compared unclassified sequences with a **reference genome** (_Enterococcus faecalis_ NZ\_OVVR01000001.1).
    -   Due to memory constraints, only **one reference sequence** was used.
4.  **Parsimony Tree Construction:**
    -   Built a parsimony-based phylogenetic tree for deeper taxonomic insights.

## **Key Findings**

-   **BLAST results (~95% identity) suggest a close relationship with _Enterococcus faecalis_.**
-   **NJ tree analysis supports genetic similarity but does not confirm identical taxonomy.**
-   **Potential evidence of a related but distinct strain or a novel taxon.**

## **Output Files**

-   **`unclassified_sequences.fasta`** → Filtered unclassified sequences.
-   **`Unclassified_sequences_nj_tree.pdf`** → Neighbor-Joining Phylogenetic Tree.
-   **`parsimony_tree.pdf`** → Parsimony-based Phylogenetic Tree.

## **Usage**

### **Reproducing the Analysis**

1.  Ensure R and required packages (`tidyverse`, `Biostrings`, `DECIPHER`, `ape`, `phangorn`) are installed.
2.  Run the R script to execute all steps.
3.  Upload `unclassified_sequences.fasta` to **NCBI BLASTn** for verification.


# **Task 2: Metagenomic Species Identification from Oxford Nanopore Data**

## **Overview**

This task identifies microbial species from an **Oxford Nanopore** metagenomic dataset (POD5 format) using **Guppy basecalling**, **Kraken2 classification**, and **Krona visualization**.

## **Workflow Summary**

1.  **Conversion:**
    -   Used `pod5 convert to_fast5` to convert **POD5** files to **FAST5** (since Guppy does not accept POD5 directly).
2.  **Basecalling:**
    -   Applied **Guppy basecaller** on FAST5 files to generate FASTQ reads.
    -   Extracted **high-quality reads** from the `pass` folder.
3.  **Taxonomic Classification:**
    -   Used **Kraken2** with a **mini database** (due to memory constraints) in **Google Colab** to classify reads.
4.  **Visualization:**
    -   Generated a **Krona chart** for taxonomic distribution analysis.

## **Results (Krona Chart Interpretation)**

-   **"No hits" (25%)** → Unclassified reads.
-   **"Pseudomonadota" (25%)** → Reads assigned to this bacterial phylum.
-   **"Methylophilales methylotrophus" (25%)** → Reads classified under this species.
-   **"Other Root" (25%)** → Miscellaneous classifications not explicitly labeled.

## **Output Files**

-   **`kraken2_report.txt`** → Kraken2 classification results.
-   **`krona_chart.html`** → Interactive Krona visualization.

## **Usage**

1.  Convert POD5 to FAST5 using.
2.  Run Guppy basecaller on FAST5 file.
3.  Classify reads with Kraken2 and visualize results using Krona

