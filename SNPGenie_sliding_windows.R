#! /usr/bin/env Rscript --slave --no-restore --no-save

############################################################################################################
# SNPGenie process codon results to bootstrap

# Mac
suppressMessages(library(package = readr))
suppressMessages(library(package = stringr))
suppressMessages(library(package = dplyr))
suppressMessages(library(package = boot))

# At command line, call something like:
# SNPGenie_sliding_windows.R codon_results_oneProduct.tsv N S 40 1 1000 6 NONE 6 > SNPGenie_sliding_windows_oneProduct.out

############################################################################################################
############################################################################################################
### GATHER GLOBAL VARIABLES FROM COMMAND LINE
ARGV <- commandArgs(trailingOnly = T)

kill_script <- FALSE

if(! (length(ARGV) >= 5)) {
  cat("\n### Command line Error 1.\n")
  kill_script <- TRUE
} else if(! str_detect(string = ARGV[1], pattern = "\\w")) {
  cat("\n### Command line Error 2.\n")
  kill_script <- TRUE
} else if(! (ARGV[2] == "N" || ARGV[2] == "T")) {
  cat("\n### Command line Error 3.\n")
  kill_script <- TRUE
} else if(! (ARGV[3] == "S" || ARGV[3] == "N")) {
  cat("\n### Command line Error 4.\n")
  kill_script <- TRUE
} else if(! (str_detect(string = ARGV[4], pattern = "\\d") && as.integer(ARGV[4]) >= 2)) {
  cat("\n### Command line Error 5.\n")
  kill_script <- TRUE
} else if(! (str_detect(string = ARGV[5], pattern = "\\d") && as.integer(ARGV[5]) >= 1)) {
  cat("\n### Command line Error 6.\n")
  kill_script <- TRUE
} # could add more conditions later

if(kill_script) {
  cat("\n\n### WARNING: there must be 5-10 command line arguments, in this order:\n")
  cat("    (1) CODON RESULTS FILE\n")
  cat("    (2) NUMERATOR SITE TYPE (\"N\" or \"T\"; e.g., \"N\" for dN)\n")
  cat("    (3) DENOMINATOR SITE TYPE (\"S\" or \"N\"; e.g., \"S\" for dS)\n")
  cat("    (4) SLIDING WINDOW SIZE (CODONS; ≥2; ≥10 recommended)\n")
  cat("    (5) SLIDING WINDOW STEP SIZE (CODONS; ≥1)\n")
  cat("    (6) NUMBER OF BOOTSTRAP REPLICATES PER WINDOW (OPTIONAL; ≥2; DEFAULT=1000)\n")
  cat("    (7) MINIMUM NUMBER OF DEFINED CODONS PER CODON POSITION (OPTIONAL; ≥2; DEFAULT=6)\n")
  cat("    (8) MULTIPLE HITS CORRECTION (OPTIONAL; \"NONE\" or \"JC\"; DEFAULT=NONE)\n")
  cat("    (9) NUMBER OF CPUS (OPTIONAL; ≥1; DEFAULT=1)\n")
  cat("    (10) STRING TO PREPEND TO OUTPUT LINES (OPTIONAL; DEFAULT=\"\")\n")
  cat("Contact author if this help is insufficient.\n\n")
  quit(save = 'no', status = 1, runLast = TRUE)
} #else {
#  cat("\n\n### RUNNING! ###\n")
#}

CODON_RESULTS_FILE <- as.character(ARGV[1])
NUMERATOR <- as.character(ARGV[2])
DENOMINATOR <- as.character(ARGV[3])
WINDOW_SIZE <- as.integer(ARGV[4]) # <- 10
STEP_SIZE <- as.integer(ARGV[5]) # <- 1
NBOOTSTRAPS <- as.integer(ARGV[6]) # 1000
MIN_DEFINED_CODONS <- as.integer(ARGV[7]) # 0, 6, or 8 # only used by SNPGenie
CORRECTION <- 'NONE' # "JC"
NCPUS <- 1 # 6
PREPEND_TO_OUTPUT <- ''

# Produce some helpful warning messages
if(! (NBOOTSTRAPS >= 1 && NBOOTSTRAPS <= 1000000)) {
  cat("\n### WARNING: NBOOTSTRAPS must be in the range [1,1000000]. Using: 1000.\n")
  NBOOTSTRAPS <- 1000
}

if(! (MIN_DEFINED_CODONS >= 2)) {
  cat("\n### WARNING: MIN_DEFINED_CODONS must be ≥2. Using: 6.\n")
  MIN_DEFINED_CODONS <- 6
}

if(! is.na(ARGV[8]) && ARGV[8] == "JC") {
  CORRECTION <- as.character(ARGV[8])
} else if (is.na(ARGV[8]) || ARGV[8] != "NONE") {
  cat("\n### WARNING: unrecognized CORRECTION supplied. Using: \"NONE\".\n")
}

if(! is.na(ARGV[9]) && str_detect(string = ARGV[9], pattern = "\\d") && ! str_detect(string = ARGV[9], pattern = "[a-zA-Z]")) {
  NCPUS <- as.integer(ARGV[9]) # <- 4, 8, or 60
}

if(! is.na(ARGV[10])) {
  PREPEND_TO_OUTPUT <- as.character(ARGV[10])
}

# EXAMPLES
#CODON_RESULTS_FILE <- "/Users/cwnelson88/Desktop/SCIENCE/SARS-CoV-2-South-Africa/intrahost_sliding_windows_size40/codon_results_byProductCodon_E_WINDOWS_dNdS.tsv"
#CODON_RESULTS_FILE <- "/Users/cwnelson88/Desktop/SCIENCE/SNPGenie/Dan/between_group_codon_results.txt" # between-group example failure
#CODON_RESULTS_FILE <- "/Users/cwnelson88/Desktop/SCIENCE/SNPGenie/Dan/between_group_codon_results_oneComp_oneProduct.txt" # between-group example success
#NUMERATOR <- "N"
#DENOMINATOR <- "S"
#WINDOW_SIZE <- 10
#STEP_SIZE <- 1
#MIN_DEFINED_CODONS <- 90 # 100 # 2
#NBOOTSTRAPS <- 100
#CORRECTION <- 'NONE'
#NCPUS <- 6
#PREPEND_TO_OUTPUT <- ''

### Read in the file
suppressMessages(codon_data <- read_tsv(file = CODON_RESULTS_FILE))

# If it's a classic snpgenie.pl output, add codon_num and num_defined_seqs columns
if('N_diffs_vs_ref' %in% names(codon_data)) {
  min_site <- min(codon_data$site)
  
  # number codons and make sure they're all evenly divisible
  if(any((codon_data$site - min_site + 3) %% 3 > 0)) {
    cat("\n###############################################################################\n")
    cat("    ERROR: the input seems to contain codons for more than one gene/product.\n")
    cat("    If this is not the case, your gene may contain multiple exons: contact author for an update.\n")
    quit(save = 'no', status = 1, runLast = TRUE)
  } else {
    codon_data$codon_num <- (codon_data$site - min_site + 3) / 3 # number the codons
  }
  
  # assume a minimum read depth of 100
  codon_data$num_defined_seqs <- 100
  cat("\n###############################################################################\n")
  cat("ALERT: bootstrapping assumes reliable estimates of N and S differences for all codons.\n")
  cat("    For snpgenie.pl results, this means there should be sufficient coverage at each site,\n")
  cat("    and variants should be filtered for false positives.\nIt is incumbent on the user to check these issues.\n")
} 

# If it's a between-group output
if('num_defined_codons_g1' %in% names(codon_data)) {
  cat("\n###############################################################################\n")
  cat("BETWEEN-GROUP FILE DETECTED. SNPGenie will require both groups to have at least\n")
  cat("the user-specified number of MIN_DEFINED_CODONS:", MIN_DEFINED_CODONS, ".\n")
  
  # Calculate minimum number defined codons for the two groups being compared
  min_num_defined_codons <- apply(data.frame(num_defined_codons_g1 = codon_data$num_defined_codons_g1, num_defined_codons_g2 = codon_data$num_defined_codons_g2),
                                  1, min) # 1 = row-wise
  codon_data$num_defined_seqs <- min_num_defined_codons
  
  # Reformat the codon_num values for consistency
  codon_data$codon_num <- as.integer(str_replace(string = codon_data$codon_num, pattern = "codon_", replacement = ""))
}

# Check for unique codon set
if(length(unique(codon_data$product)) != 1 || length(codon_data$codon_num) != length(unique(codon_data$codon_num))) {
  cat("\n\n### WARNING: there must be results for ONLY ONE GENE PRODUCT PER FILE, i.e., a unique set of codons per frame.\n")
  cat("### Additionally, for between-group comparisons, there must be results for ONLY ONE GROUP-TO-GROUP COMPARISON.\n")
  cat("### TERMINATED.\n\n")
  quit(save = 'no', status = 1, runLast = TRUE)
}

### PRINT ANALYSIS LOG:
cat("\n###############################################################################\n")
cat("SNPGenie sliding window analysis LOG\n\n")
cat("PARAMETERS:\n")
cat(paste0("(1) CODON_RESULTS_FILE=", CODON_RESULTS_FILE, "\n"))
cat(paste0("(2) NUMERATOR=", NUMERATOR, "\n"))
cat(paste0("(3) DENOMINATOR=", DENOMINATOR, "\n"))
cat(paste0("(4) WINDOW_SIZE=", WINDOW_SIZE, "\n"))
cat(paste0("(5) STEP_SIZE=", STEP_SIZE, "\n"))
cat(paste0("(6) NBOOTSTRAPS=", NBOOTSTRAPS, "\n"))
cat(paste0("(7) MIN_DEFINED_CODONS=", MIN_DEFINED_CODONS, "\n"))
cat(paste0("(8) CORRECTION=", CORRECTION, "\n"))
cat(paste0("(9) NCPUS=", NCPUS, "\n"))
cat(paste0("(10) PREPEND_TO_OUTPUT=", PREPEND_TO_OUTPUT, "\n\n"))


############################################################################################################
############################################################################################################
### BOOTSTRAP FUNCTIONS
############################################################################################################
############################################################################################################

############################################################################################################
# BOOTSTRAP FUNCTION (dN - dS) for CODON UNIT, NO CORRECTION
dNdS_diff_boot_fun <- function(codon_results, numerator, denominator, num_replicates, num_cpus) {
  
  # Function for dN
  dN_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    return(dN)
  }
  
  # Function for dN
  dS_function <- function(D, indices) {
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    return(dS)
  }
  
  # Function for dN - dS
  dN_m_dS_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dN_m_dS <- dN - dS
    return(dN_m_dS)
  }
  
  # Function for dN/dS
  dN_over_dS_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dN_over_dS <- dN / dS
    return(dN_over_dS)
  }
  
  # CREATE FUNCTION FOR dN/dS TO CALCULATE ITS SE
  (dN <- sum(as.vector(codon_results[ , paste0(numerator, "_diffs")]), na.rm = TRUE) / sum(as.vector(codon_results[ , paste0(numerator, "_sites")]), na.rm = TRUE))
  (dS <- sum(as.vector(codon_results[ , paste0(denominator, "_diffs")]), na.rm = TRUE) / sum(as.vector(codon_results[ , paste0(denominator, "_sites")]), na.rm = TRUE))
  (dNdS <- dN / dS)
  
  # Run the BOOTSTRAPS
  # boot dN
  (boot_dN <- boot(data = codon_results, R = num_replicates, statistic = dN_function, parallel = 'multicore', ncpus = num_cpus))
  (dN <- boot_dN$t0)
  (boot_dN_SE <- sd(boot_dN$t))
  
  # boot dS
  (boot_dS <- boot(data = codon_results, R = num_replicates, statistic = dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dS <- boot_dS$t0)
  (boot_dS_SE <- sd(boot_dS$t))
  
  # boot dN - dS
  (boot_dN_m_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_m_dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dN_m_dS <- boot_dN_m_dS$t0)
  (boot_dN_m_dS_SE <- sd(boot_dN_m_dS$t))
  (boot_dN_m_dS_Z <- dN_m_dS / boot_dN_m_dS_SE)
  (boot_dN_m_dS_P <- 2 * pnorm(-abs(boot_dN_m_dS_Z)))
  
  # boot dN/dS
  (boot_dN_over_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_over_dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dN_over_dS <- boot_dN_over_dS$t0)
  (boot_dN_over_dS_SE <- sd(boot_dN_over_dS$t))
  (boot_dN_over_dS_Z <- dN_over_dS / boot_dN_over_dS_SE)
  (boot_dN_over_dS_P <- 2 * pnorm(-abs(boot_dN_over_dS_Z)))
  
  ### NEW: ASL (acheived significance level)
  boot_dN_gt_dS_count <- sum(boot_dN_m_dS$t > 0)
  boot_dN_eq_dS_count <- sum(boot_dN_m_dS$t == 0)
  boot_dN_lt_dS_count <- sum(boot_dN_m_dS$t < 0)
  ASL_dN_gt_dS_P <- boot_dN_lt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
  ASL_dN_lt_dS_P <- boot_dN_gt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
  
  return(paste(num_replicates, dN, dS, dNdS, dN_m_dS, boot_dN_SE, boot_dS_SE, boot_dN_over_dS_SE, boot_dN_over_dS_P, 
               boot_dN_m_dS_SE, boot_dN_m_dS_P, 
               boot_dN_gt_dS_count, boot_dN_eq_dS_count, boot_dN_lt_dS_count, ASL_dN_gt_dS_P, ASL_dN_lt_dS_P,
               sep = "\t"))
}


############################################################################################################
# BOOTSTRAP FUNCTION (dN - dS) for CODON UNIT, ***JUKES-CANTOR CORRECTION***
dNdS_diff_boot_fun_JC <- function(codon_results, numerator, denominator, num_replicates, num_cpus) {
  
  # Function for dN
  dN_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dN <- -3/4 * log(1 - (4/3) * dN)
    return(dN)
  }
  
  # Function for dN
  dS_function <- function(D, indices) {
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dS <- -3/4 * log(1 - (4/3) * dS)
    return(dS)
  }
  
  # Function for dN - dS
  dN_m_dS_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dN <- -3/4 * log(1 - (4/3) * dN)
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dS <- -3/4 * log(1 - (4/3) * dS)
    dN_m_dS <- dN - dS
    return(dN_m_dS)
  }
  
  # Function for dN/dS
  dN_over_dS_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dN <- -3/4 * log(1 - (4/3) * dN)
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dS <- -3/4 * log(1 - (4/3) * dS)
    dN_over_dS <- dN / dS
    return(dN_over_dS)
  }
  
  # CREATE FUNCTION FOR dN/dS TO CALCULATE ITS SE
  dN <- sum(as.vector(codon_results[ , paste0(numerator, "_diffs")]), na.rm = TRUE) / sum(as.vector(codon_results[ , paste0(numerator, "_sites")]), na.rm = TRUE)
  dN <- -3/4 * log(1 - (4/3) * dN)
  dS <- sum(as.vector(codon_results[ , paste0(denominator, "_diffs")]), na.rm = TRUE) / sum(as.vector(codon_results[ , paste0(denominator, "_sites")]), na.rm = TRUE)
  dS <- -3/4 * log(1 - (4/3) * dS)
  dNdS <- dN / dS
  
  # Run the BOOTSTRAPS
  # boot dN
  (boot_dN <- boot(data = codon_results, R = num_replicates, statistic = dN_function, parallel = 'multicore', ncpus = num_cpus))
  (dN <- boot_dN$t0)
  (boot_dN_SE <- sd(boot_dN$t))
  
  # boot dS
  (boot_dS <- boot(data = codon_results, R = num_replicates, statistic = dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dS <- boot_dS$t0)
  (boot_dS_SE <- sd(boot_dS$t))
  
  # boot dN - dS
  (boot_dN_m_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_m_dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dN_m_dS <- boot_dN_m_dS$t0)
  (boot_dN_m_dS_SE <- sd(boot_dN_m_dS$t))
  (boot_dN_m_dS_Z <- dN_m_dS / boot_dN_m_dS_SE)
  (boot_dN_m_dS_P <- 2 * pnorm(-abs(boot_dN_m_dS_Z)))
  
  # boot dN/dS
  (boot_dN_over_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_over_dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dN_over_dS <- boot_dN_over_dS$t0)
  (boot_dN_over_dS_SE <- sd(boot_dN_over_dS$t))
  (boot_dN_over_dS_Z <- dN_over_dS / boot_dN_over_dS_SE)
  (boot_dN_over_dS_P <- 2 * pnorm(-abs(boot_dN_over_dS_Z)))
  
  ### NEW: ASL (acheived significance level)
  boot_dN_gt_dS_count <- sum(boot_dN_m_dS$t > 0)
  boot_dN_eq_dS_count <- sum(boot_dN_m_dS$t == 0)
  boot_dN_lt_dS_count <- sum(boot_dN_m_dS$t < 0)
  ASL_dN_gt_dS_P <- boot_dN_lt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
  ASL_dN_lt_dS_P <- boot_dN_gt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
  
  return(paste(num_replicates, dN, dS, dNdS, dN_m_dS, boot_dN_SE, boot_dS_SE, boot_dN_over_dS_SE, boot_dN_over_dS_P, 
               boot_dN_m_dS_SE, boot_dN_m_dS_P, 
               boot_dN_gt_dS_count, boot_dN_eq_dS_count, boot_dN_lt_dS_count, ASL_dN_gt_dS_P, ASL_dN_lt_dS_P,
               sep = "\t"))
}

############################################################################################################
############################################################################################################


############################################################################################################
### PERFORM SLIDING WINDOWS WITH BOOTSTRAPS

### ADD COLUMNS TO DATA
RATIO_NAME <- paste0('d', NUMERATOR, 'd', DENOMINATOR)
codon_data$sw_ratio <- RATIO_NAME
codon_data$sw_start <- as.double(NA)
codon_data$sw_center <- as.double(NA)
codon_data$sw_end <- as.double(NA)
codon_data$sw_num_replicates <- as.integer(NA)
codon_data$sw_N_diffs <- as.double(NA)
codon_data$sw_S_diffs <- as.double(NA)
codon_data$sw_N_sites <- as.double(NA)
codon_data$sw_S_sites <- as.double(NA)
codon_data$sw_dN <- as.double(NA)
codon_data$sw_dS <- as.double(NA)
codon_data$sw_dNdS <- as.double(NA)
codon_data$sw_dN_m_dS <- as.double(NA)
codon_data$sw_boot_dN_SE <- as.double(NA)
codon_data$sw_boot_dS_SE <- as.double(NA)
codon_data$sw_boot_dN_over_dS_SE <- as.double(NA)
codon_data$sw_boot_dN_over_dS_P <- as.double(NA)
codon_data$sw_boot_dN_m_dS_SE <- as.double(NA)
codon_data$sw_boot_dN_m_dS_P <- as.double(NA)
codon_data$sw_boot_dN_gt_dS_count <- as.integer(NA)
codon_data$sw_boot_dN_eq_dS_count <- as.integer(NA)
codon_data$sw_boot_dN_lt_dS_count <- as.integer(NA)
codon_data$sw_ASL_dN_gt_dS_P <- as.double(NA)
codon_data$sw_ASL_dN_lt_dS_P <- as.double(NA)
codon_data$sw_ASL_dNdS_P <- as.double(NA)


### PERFORM SLIDING WINDOW
cat("Performing sliding window analysis...\n\n")

codon_data_filtered <- filter(codon_data, num_defined_seqs >= MIN_DEFINED_CODONS)
codon_data_filtered <- dplyr::arrange(codon_data_filtered, codon_num)

if(nrow(codon_data_filtered) >= WINDOW_SIZE) {
  for(i in (min(codon_data_filtered$codon_num)):(max(codon_data_filtered$codon_num) - WINDOW_SIZE + 1)) {
    #for(i in 1:(nrow(codon_data_filtered) - WINDOW_SIZE + 1)) { # each window starting at row 1
    #i <- 1 # 111
    #cat(i, "")
    
    # Extract window; analyze
    window_codon_data <- codon_data_filtered[i:(i + WINDOW_SIZE - 1), ]
    lowest_codon_num <- min(window_codon_data$codon_num)
    highest_codon_num <- max(window_codon_data$codon_num)
    
    if(nrow(window_codon_data) > 1) {
      
      # WARNING if missing codon numbers
      if(nrow(window_codon_data) != WINDOW_SIZE) {
        cat("### WARNING: codon", i, "is missing.")
        cat("### If this was unexpected or is common, there may be problems with your data.")
      }
      
      ### BOOTSTRAP EACH RATIO FOR THIS WINDOW. Columns:
      #[1] num_replicates
      #[2] dN
      #[3] dS
      #[4] dNdS
      #[5] dN_m_dS
      #[6] boot_dN_SE
      #[7] boot_dS_SE
      #[8] boot_dN_over_dS_SE
      #[9] boot_dN_over_dS_P
      #[10] boot_dN_m_dS_SE
      #[11] boot_dN_m_dS_P
      #[12] boot_dN_gt_dS_count
      #[13] boot_dN_eq_dS_count
      #[14] boot_dN_lt_dS_count
      #[15] ASL_dN_gt_dS_P
      #[16] ASL_dN_lt_dS_P
      
      # dN/dS bootstrap with appropriate correction
      boot_dNdS <- NA
      if(CORRECTION == "JC") {
        boot_dNdS <- dNdS_diff_boot_fun_JC(window_codon_data, NUMERATOR, DENOMINATOR, NBOOTSTRAPS, NCPUS)
      } else {
        boot_dNdS <- dNdS_diff_boot_fun(window_codon_data, NUMERATOR, DENOMINATOR, NBOOTSTRAPS, NCPUS)
      }
      boot_dNdS <- str_split(string = boot_dNdS, pattern = '\t')[[1]]
      
      # bootstrap results
      num_replicates <- as.integer(boot_dNdS[1])
      dN <- as.numeric(boot_dNdS[2])
      dS <- as.numeric(boot_dNdS[3])
      dNdS <- as.numeric(boot_dNdS[4])
      dN_m_dS <- as.numeric(boot_dNdS[5])
      suppressWarnings(boot_dN_SE <- as.numeric(boot_dNdS[6]))
      suppressWarnings(boot_dS_SE <- as.numeric(boot_dNdS[7]))
      suppressWarnings(boot_dN_over_dS_SE <- as.numeric(boot_dNdS[8]))
      suppressWarnings(boot_dN_over_dS_P <- as.numeric(boot_dNdS[9]))
      suppressWarnings(boot_dN_m_dS_SE <- as.numeric(boot_dNdS[10]))
      suppressWarnings(boot_dN_m_dS_P <- as.numeric(boot_dNdS[11]))
      suppressWarnings(boot_dN_gt_dS_count <- as.integer(boot_dNdS[12]))
      suppressWarnings(boot_dN_eq_dS_count <- as.integer(boot_dNdS[13]))
      suppressWarnings(boot_dN_lt_dS_count <- as.integer(boot_dNdS[14]))
      suppressWarnings(ASL_dN_gt_dS_P <- as.numeric(boot_dNdS[15]))
      suppressWarnings(ASL_dN_lt_dS_P <- as.numeric(boot_dNdS[16]))
      ASL_dNdS_P <- 1
      
      if(! is.na(ASL_dN_gt_dS_P) && ! is.na(ASL_dN_lt_dS_P) && ASL_dN_gt_dS_P < ASL_dN_lt_dS_P) {
        ASL_dNdS_P <- 2 * ASL_dN_gt_dS_P
      } else if(! is.na(ASL_dN_gt_dS_P) && ! is.na(ASL_dN_lt_dS_P)) {
        ASL_dNdS_P <- 2 * ASL_dN_lt_dS_P
      }
      
      if(ASL_dNdS_P == 0) {
        ASL_dNdS_P <- 1 / NBOOTSTRAPS
      }
      
      # Add to table
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_start <- lowest_codon_num
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_center <- (lowest_codon_num + highest_codon_num) / 2
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_end <- highest_codon_num
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_num_replicates <- num_replicates
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_N_diffs <- sum(unname(unlist(window_codon_data[ , paste0(NUMERATOR, '_diffs')])), na.rm = TRUE)
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_N_sites <- sum(unname(unlist(window_codon_data[ , paste0(DENOMINATOR, '_diffs')])), na.rm = TRUE)
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_S_diffs <- sum(unname(unlist(window_codon_data[ , paste0(NUMERATOR, '_sites')])), na.rm = TRUE)
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_S_sites <- sum(unname(unlist(window_codon_data[ , paste0(DENOMINATOR, '_sites')])), na.rm = TRUE)
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_dN <- dN
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_dS <- dS
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_dNdS <- dNdS
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_dN_m_dS <- dN_m_dS
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dN_SE <- boot_dN_SE
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dS_SE <- boot_dS_SE
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dN_over_dS_SE <- boot_dN_over_dS_SE
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dN_over_dS_P <- boot_dN_over_dS_P
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dN_m_dS_SE <- boot_dN_m_dS_SE
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dN_m_dS_P <- boot_dN_m_dS_P
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dN_gt_dS_count <- boot_dN_gt_dS_count
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dN_eq_dS_count <- boot_dN_eq_dS_count
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_boot_dN_lt_dS_count <- boot_dN_lt_dS_count
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_ASL_dN_gt_dS_P <- ASL_dN_gt_dS_P
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_ASL_dN_lt_dS_P <- ASL_dN_lt_dS_P
      codon_data[codon_data$codon_num == lowest_codon_num, ]$sw_ASL_dNdS_P <- ASL_dNdS_P
      
    } # else leave it as NA
  } # end last window
} else if(nrow(codon_data) < WINDOW_SIZE) {
  cat("\n\n### WARNING: the number of codons with num_defined_seqs>", MIN_DEFINED_CODONS , " is less than the WINDOW_SIZE=", WINDOW_SIZE, ".\n", sep="")
  cat("### Thus, a sliding window using these parameters is not possible.\n")
  cat("### TERMINATED.\n\n")
  quit(save = 'no', status = 1, runLast = TRUE)
  
} # end if at least WINDOW_SIZE

### SAVE RESULTS
OUTFILE <- CODON_RESULTS_FILE
OUTFILE <- str_replace(string = OUTFILE, pattern = ".txt", replacement = "")
OUTFILE <- str_replace(string = OUTFILE, pattern = ".tsv", replacement = "")
OUTFILE <- paste0(OUTFILE, "_WINDOWS_", RATIO_NAME, ".tsv")

cat(paste0("Writing output to: ", OUTFILE, "\n\n"))

write_tsv(codon_data, OUTFILE)

cat("DONE\n\n")

