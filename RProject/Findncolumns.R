base.freq(dna_bin[,1],all =TRUE)
f1 <- function(col1){out1 <- base.freq(col1, all = TRUE)
  out1 <- out1["n"]>0
  return(out1)
}

keep_col <- NULL
for(i in 1:ncol(dna_bin)){keep_col <- c(keep_col, f1(dna_bin[,i]))}
reduc1 <- dna_bin[,!keep_col]


# Calculate the Hamming distance matrix
hamming_dist <- dist.hamming(reduc1)
# Print the distance matrix
print(mean(hamming_dist))



# Calculate the pairwise differences
pairwise_diff <- dist.dna(reduc1, model = "n", pairwise.deletion = FALSE)
# Print the pairwise differences
print(mean(pairwise_diff))