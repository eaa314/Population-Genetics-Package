library(pinfsc50)
#Locate (get path to) file
system.file("extdata","pinf_sc50.vcf.gz", package = "pinfsc50")

library(vcfR)
library(ape)
# Set the path to your VCF file
vcf_file <- "/home/eaa31/R/x86_64-pc-linux-gnu-library/4.3/pinfsc50/extdata/pinf_sc50.vcf.gz"
# Read the VCF file
vcf <- read.vcfR(vcf_file)
# Convert VCF to DNAbin format
dna_bin <- vcfR2DNAbin(vcf)
# Save the DNAbin object to a file
#write.dna(dna_bin, file = "output.dna", format = "fasta")

#all(nrow(dna_bin)*base.freq(dna_bin[,1])!=1)



Achaz_correction <- function(dna_bin) {

  # Initialize an empty list to store selected columns
  #selected_columns <- list()

  selected_columns <- NULL

  # Iterate over each column in the DNAbin object
  for (col_index in 1:ncol(dna_bin)) {
    # Check the condition for the column
    if (all(nrow(dna_bin) * base.freq(dna_bin[, col_index], all = TRUE)[c("a","c","t","g")] != 1)) {
      # If condition is true, store the column in the list
      #selected_columns[[length(selected_columns) + 1]] <- dna_bin[, col_index]
      selected_columns <- c(selected_columns, TRUE)
    } else{selected_columns <- c(selected_columns, FALSE)}
  }
  str(selected_columns)

  # Combine selected columns into a new DNAbin object
  #new_dna_bin <- do.call(cbind, selected_columns)
  new_dna_bin <- dna_bin[,selected_columns]
  str(new_dna_bin)

  # Print the new DNAbin object
  print(new_dna_bin)
  View(new_dna_bin)

  #Read the DNAbin file
  #dna_bin <- read.dna("/home/eaa31/R_Analyses/output.dna", format = "fasta")

#  # Calculate the Hamming distance matrix
#  hamming_dist <- dist.hamming(dna_bin)
#  # Print the distance matrix
#  print(mean(hamming_dist))

  # Calculate the pairwise differences
  pairwise_diff <- dist.dna(new_dna_bin, model = "n", pairwise.deletion = FALSE)
  # Print the pairwise differences
  #print(mean(pairwise_diff))

#  # Convert DNAbin to matrix representation
#  dna_matrix <- as.matrix(dna_bin)
#  # Count the number of segregating sites
#  segregating_sites <- sum(apply(dna_matrix, 2, function(x) length(intersect(unique(as.character(x)), c("a","c","t","g"))) > 1))
#  # Print the number of segregating sites
#  print(segregating_sites)

  # Calculate the number of segregating sites
  segregating_sites <- seg.sites(new_dna_bin)
  
  # Print the result
  #print("Number of segregating sites:")
  #print(length(segregating_sites))

  #sw(dna_bin, 100, 100)

  #p = pairwise differences (pi)
  #s = segregating sites
  #t1 = theta pi-n1
  #e1 = estimated t1
  #t2 = theta s-n1
  #e2 = estimated t2

  p= mean(pairwise_diff)
  n= nrow(dna_bin)
  s= length(segregating_sites)

  #t1 <- p * (n - 1) / (n - 3)
  #e1 <- t1 * (n - 3) / (n - 1)
  #t2 <- (s - n) / (an - n / (n - 1))
  #e2 <- t2 * (an - (n / (n - 1)))
  
  #Formula B1
  an <- sum(1/(1:(n-1)))
  bn <- sum(1/(1:(n-1))^2)
  
  #Gamma n and sigma n, within formula B22
  g <- an - (n / (n - 1))
  si <- bn - ((2*n - 1) / (n - 1)^2)
  
  #Theta and Theta^2, fourmulae B25 and B26
  t <- s / g
  t2 <- (s * (s - 1)) / (si + g^2)
  
  #f*, formula 22
  f <- (n - 3) / (an * (n - 1) - n)

  #Formulae B20 and B21
  alphan <- f^2 * (an - (n) / (n - 1)) + f * (an * (4*(n + 1)/(n - 1)^2 - 2*(n + 3)/(n - 1)) - an * (8*(n + 1)/(n*(n - 1)^2)) + (n^2 + n + 60)/(3*n*(n - 1)))
  betan <- f^2 * (bn - (2*n - 1)/(n - 1)^2) + f * (bn * (8/(n - 1)) - an * (4/(n*(n - 1))) - (n^3 + 12*n^2 - 35*n + 18)/(n*(n - 1)^2)) - bn * (16/(n*(n - 1)) + an * (8/(n^2*(n - 1))) + 2*(n^4 + 110*n^2 - 255*n + 126)/(9*n^2*(n - 1)^2))

  #Formula B19
  variance <- alphan*t - betan*t2
  
  #Theta estimators s-n1 and pi-n1, formulae 13 and 15
  #e1 <- s / (an - n / (n - 1))
  #e2 <- p * (n - 1) / (n - 3) 

  y <- (p - f * s) / sqrt(variance)
  return(y)
}

Achaz_correction(new_dna_bin)
