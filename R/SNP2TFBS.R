findMotifs <- function(pwmList, input_high, input_low, background_1, background_2, seq_width, percentage = "85%", output_file = F, test_run = F) {
  
  ### This function will take sequences of the same regions that have high or low signal from ATAC/ChIP-seq between different alleles and calculate significant enrichment of motifs compared to background regions with low ATAC/ChIP-seq signal.
  ### Assumptions: the regions in the input high and low MUST be the same. The regions in the background_1 and background_2 MUST be the same. 
  ###### The number and sizes of high and low regions should be the same (set by seq_width). If there are more regions or they are bigger than seq_width, will try to reduce to the same number.
  
  ################## Defines function to import fasta files
  
  importFasta <- function(x) {
    inputFasta <- Biostrings::DNAStringSet()
    
    for (file in x) {
      fasta_file <- Biostrings::readDNAStringSet(file)
      inputFasta <- c(inputFasta, fasta_file)
    } 
    
    inputFasta <- inputFasta[order(names(inputFasta)),]
    #Order sequences by their region coordinates
    
    if (test_run != F & test_run %% 1 == 0) {
      
      print(paste("Testing mode on, subsetting the first", test_run, "regions from fasta!!!"))
      
      inputFasta <- inputFasta[1:test_run]
      
    }
    
    
    ###### Checks if the sizes of the input regions are equivalent to seq_width from the function
    ###### If more than one region sizes, it tries to reduce down to seq_width. If not possible, throws error
    ###### If only one region size, but larger than seq_width, resizes. Otherwise, error.
    
    
    if (length(unique(GenomicRanges::width(inputFasta))) > 1) {
      if (min(GenomicRanges::width(inputFasta)) >= seq_width) {
        print(paste("Number of nucleotides is different across", deparse(substitute(x)), "sequences! Resizing to", seq_width))
        inputFasta <- Biostrings::subseq(inputFasta, start = 1, end = seq_width)
      } else {
        stop(paste("ERROR: Some regions are smaller than", seq_width, "bp! Can't proceed."))
      }
    } else if (length(unique(GenomicRanges::width(inputFasta))) == 1) {
      if (unique(GenomicRanges::width(inputFasta)) > seq_width) {
        print(paste("Region sizes are larger than", seq_width, "bp! Downsizing to", seq_width, "bp..."))
        inputFasta <- Biostrings::subseq(inputFasta, start = 1, end = seq_width)
      } else if (unique(GenomicRanges::width(inputFasta)) == seq_width) {
        print(paste("Regions have the expected sizes, continuing..."))
      } else {
        stop(paste("ERROR: Some regions are smaller than", seq_width, "bp! Can't proceed."))
      }
    } else{
      stop(paste("ERROR: Can't proceed with", deparse(substitute(x))))
    }
    
    
    return(inputFasta)
  }
  
  ################## Imports input files with high signal
  
  inputFasta_high <- importFasta(input_high)
  
  ################## Imports input files with low signal
  
  inputFasta_low <- importFasta(input_low)
  
  ########## Checks if the number of input sequences in the high and low groups are the same
  
  if (length(names(inputFasta_high)) != length(names(inputFasta_low))) {
    stop("Number of sequences in input high are different than input low")
  }
  
#### Function to calculate motifs from .fasta file(s)
  
  calc_motifs <- function(fasta, name_col) {
    
    print(paste0("Searching for motifs within ", name_col, "..."))
    
    motifs <- searchSeq(pwmList, fasta, min.score=percentage, strand="*") %>%
      as.data.frame(writeGFF3(., scoreType = "relative")) %>%
      group_by(ID, TF) %>%
      summarise(n_motif = n(), .groups = "drop_last") %>%
      ungroup() %>%
      arrange(desc(n_motif))
    
    names(motifs)[names(motifs) == "n_motif"] <- name_col
    
    return(motifs)
    
  }
  
  ################## Calculate motifs present in input files with high signal
  
  input_motifs_high <- calc_motifs(fasta = inputFasta_high, name_col = "input_high")
  
  #Extracts the number of times a TF shows up in the results for the 'High' file. 
  
  print(paste("Found a total of", n_distinct(input_motifs_high$TF), "different TF motifs within input sequences with high signal."))
  
  ################## Calculate motifs present in input files with low signal
  
  input_motifs_low <- calc_motifs(fasta = inputFasta_low, name_col = "input_low")
  
  #Extracts the number of times a TF shows up in the results for the 'Low' file. 
  
  print(paste("Found a total of", n_distinct(input_motifs_low$TF), "different TF motifs within input sequences with low signal."))
  
  ################## Calculates the difference between the number of times a TF shows up in the 'high' vs. 'low' input files.  
  
  print(paste("Finding motifs enriched in high vs low input signal..."))
  
  motifs_enrich <- full_join(input_motifs_high, input_motifs_low, by = c("ID", "TF"))
  
  ################## Imports background_1 sequences
  
  backgroundFasta_1 <- importFasta(background_1)
  
  if (length(names(backgroundFasta_1)) != length(names(inputFasta_high))) {
    print(paste("Number of background_1 sequences is different than high input sequences, randomly picking", length(names(inputFasta_high)), "sequences"))
    backgroundFasta_1[sample(names(backgroundFasta_1), length(names(inputFasta_high)))]
  }
  #Checks if the number of background_1 sequences is the same as input_high sequences, otherwise randomly pick n sequences
  
  
  ################## Imports background_2 sequences
  
  backgroundFasta_2 <- importFasta(background_2)
  
  if (length(names(backgroundFasta_2)) != length(names(inputFasta_low))) {
    print(paste("Number of background_2 sequences is different than low input sequences, randomly picking", length(names(inputFasta_low)), "sequences"))
    backgroundFasta_2[sample(names(backgroundFasta_2), length(names(inputFasta_low)))]
  }
  #Checks if the number of background_2 sequences is the same as input_low sequences, otherwise randomly pick n sequences
  
  ################## Calculate motifs present within background_1
  
  background_motifs_1 <- calc_motifs(fasta = backgroundFasta_1, name_col = "background_1")
  
  #Extracts the number of times a TF shows up in the results from background_1
  
  print(paste("Found a total of", n_distinct(background_motifs_1$TF), "different TF motifs within background_1 sequences."))
  
  ################## Calculate motifs present within background_2
  
  background_motifs_2 <- calc_motifs(fasta = backgroundFasta_2, name_col = "background_2")
  
  #Extracts the number of times a TF shows up in the results from background_2
  
  print(paste("Found a total of", n_distinct(background_motifs_2$TF), "different TF motifs within background_2 sequences."))
  
  ################## Calculates the difference between the number of times a TF shows up in background_1 vs background_2.
  
  print(paste("Finding motifs enriched in background_1 vs background_2..."))
  
  background_enrich <- full_join(background_motifs_1, background_motifs_2, by = c("ID", "TF"))
  
  ############ Merge input and background tables together.
  
  input_background_enrich <- full_join(motifs_enrich, background_enrich, by = c("ID", "TF")) %>%
    mutate(across(c(input_high, input_low, background_1, background_2), ~replace_na(.x, 0)))
  
  ############ Calculate fisher's exact test
  
  enrich_pvalue <- input_background_enrich %>%
    select(ID, TF, input_high, input_low, background_1, background_2) %>%
    mutate(input_diff = input_high - input_low,
           background_diff = background_1 - background_2) %>%
    #Gets the difference between input high and low and background_1 and background_2
    filter(!if_all(contains("diff"), ~ . == 0)) %>%
    #Filters rows in which the difference is zero for both input and background
    rowwise %>%
    mutate(input_max = max(input_high, input_low),
           background_max = max(background_1, background_2)) %>%
    #Creates a column with the highest count of motifs in either high or low regions
    #The logic here is that if a motif is gained in high vs. low, it will use the total amount of times a motif is found within the high sequences
    #If a motif is lost in high vs. low, it will use the total amount of times a motif is found within the low sequences 
    mutate(p=fisher.test(matrix(c(abs(input_diff), input_max, abs(background_diff), background_max), nrow=2))$p.value) %>%
    #Fisher test of the absolute difference between high vs. low and the highest amount of motifs in either high or low
    #The logic here is use the difference between high and low over the total amount of motifs in either high or low, depending if a motif is overall gained or lost
    ungroup() %>%
    mutate(p.adj = p.adjust(.$p, method = "BH", n = length(.$p))) %>%
    rstatix::add_significance(p.col = "p.adj",
                              output.col = "p.adj.signif",
                              cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                              symbols = c("***", "**", "*", "ns"))
  
  enrich_diff <- enrich_pvalue %>%
    select(ID, TF, input_high, input_low, input_diff, background_1, background_2, background_diff, p, p.adj, p.adj.signif) %>%
    arrange(p.adj)
  
  ### Saves file if output_file is not F
  if (output_file != F) {
    write_csv(enrich_diff, file = output_file)
    print(paste("Exported results to:", output_file))
    invisible()
  } else {
    return(enrich_diff)
  }
  
}