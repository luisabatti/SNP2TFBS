# SNP2TFBS
An R script that identifies the enrichment or loss of transcription factor binding sites at single nucleotide polymorphisms (SNPs) associated with changes in ChIP-seq/ATAC-seq across multiple genomic regions.

# Usage

👉 You can find a step-by-step tutorial on how to use this script here: https://luisabatti.github.io/SNP2TFBS/

👇 Here are the things that you will need to run this script:

1. A list of Position-specific Weight Matrices (PWMs): I recommend using JASPAR2022.
2. Two sets of input .fasta files containing the regions to be compared. One file must contain the sequences that show high ChIP-seq/ATAC-seq signal. The other file must contain the sequences that show low ChIP-seq/ATAC-seq signal. They must be the ** SAME ** regions but with and without a SNP (for example, WT and mutated sequences). 
3. Two sets of background .fasta files containing regions to be used as a control. Like the input set, you must use the ** SAME ** regions, with and without a SNP. However, over here you should use regions that do NOT show any gain/loss of ChIP-seq/ATAC-seq signal. 

## Required libraries

👉 This script requires the following libraries:

1. tidyverse
2. TFBSTools
3. JASPAR2022

# Input

👉 The script takes the following parameters: 

1. pwmList: PWM list from JASPAR2022.
2. input_high: .fasta file containing input sequences with high ChIP-seq/ATAC-seq signal. 
3. input_low: .fasta file containing input sequences with low ChIP-seq/ATAC-seq signal. 
4. background_1: .fasta file containing background sequences, must contain the same number of sequences as input_high.   
5. background_2: .fasta file containing background sequences, must contain the same number of sequences as input_low.
6. seq_width: the length of each sequence within each .fasta file. 
7. percentage: the threshold for detecting a motif in TFBSTools. Default: 85%.
8. output_file: path to output .csv file. If set to False, return the results to a variable. Default: False.
9. test_run: subset the .fasta files to N regions to test the script. Default: False.

# Output

👇 Here are the different columns that you will find in the output .csv file:

1. ID: JASPAR2022 motif ID.
2. TF: Transcription factor name.
3. input_high: number of times a motif appeared within input_high sequences.
4. input_low: number of times a motif appeared within input_low sequences.
5. input_diff: the difference in the number of motifs between input_high and input_low.
6. background_1: number of times a motif appeared within background_1 sequences.
7. background_2: number of times a motif appeared within background_2 sequences.
8. background_diff: the difference in the number of motifs between background_1 and background_2.
9. p: Fisher's exact test of the absolute difference between high vs. low and the highest amount of motifs in either high or low. 
10. p.adj: FDR-adjusted p value.
11. p.adj.signif: Significance level (*** < 0.001, ** < 0.01, * < 0.05, ns = non-significant).

# Credits

Written by Dr. Luis E. Abatti