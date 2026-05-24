#' Omega Plus Statistic
#'
#' Calculates LD-Based Omega Statistic
#' @references Kim, Y. and Nielsen, R. (2004) ‘Linkage disequilibrium as a signature of selective sweeps’, Genetics, 167(3), pp. 1513–1524. doi:10.1534/genetics.103.025387. 
#' @param dna_bin Genomic dataset as a DNA bin file, define position information when calling function
#' @return Omega Statistic as numerical values in a data frame
#' @examples 
#' Omega_statistic(dna_bin)
#' @export
Omega_statistic <- function(dna_bin) {

  data_test <- dna_bin
  colnames(data_test) <- paste0("snp",1:ncol(data_test))
  #data_test <- as.matrix(data_test)
  data_test <- as.character(data_test)

  #data_test <- read.dna("test_data.fasta",format = "fasta",
  #                      as.character = TRUE)

  #Function to deal with NAs in the dataset
  recode_dna <- function(col1){
  #col1 <- as.character(col1)
  seq_names <- names(col1)
  tab1 <- table(col1)
  bestnuc <- names(which.max(tab1))
  nucset <- c("a","c","t","g")
  nucset2 <- setdiff(nucset,bestnuc)
  col1[!(col1 %in% nucset)] <- NA
  col1[col1==bestnuc] <- 1
  col1[col1 %in% nucset2] <- 0
  #apply(col1,1,as.integer)
  out1 <- as.integer(col1)
  names(out1) <- seq_names
  return(out1)
  }

  recoded_matrix <- apply(data_test,2,recode_dna)
  #str(recoded_matrix)

  #Creates a matrix of correlation coefficients
  correlation_coefficient <- cor(recoded_matrix, use = "pairwise.complete.obs")
  #correlation_coefficient[,1:10]

  #summary(as.vector(correlation_coefficient[1:1000,1:1000]))
  #View(correlation_coefficient[31:50,31:50])
  #View(data_test[,"snp35"])

  #Preparing datasets and variables to carry out Omega Plus calculations
  correlation_coefficient <- correlation_coefficient^2
  diag(correlation_coefficient) <- 0
  snp_names <- colnames(recoded_matrix)
  S <- ncol(correlation_coefficient)
  omega_values <- data.frame(l = numeric(), omega = numeric())

  gc()

#For loop that finds omega values using each matrix entry as a cutoff and compiles them into a single table
  for (l in 1:S) {
  #for (l in 1:10) {
    #l <- 6
    L <- snp_names[1:l]
    R <- snp_names[-c(1:l)]

    # Preparing values to use in Omega formula
    L_values <- correlation_coefficient[L,L]
    #L_values <- L_values[lower.tri(L_values)]
    #L_values <- (sum(L_values, na.rm = TRUE) - sum(diag(L_values), na.rm = TRUE)) / 2
    #Lr2 <- sum(L_values,na.rm = TRUE)
    Lr2 <- sum(L_values, na.rm = TRUE) /2
    Lnna <- sum(!is.na(L_values))/2

    R_values<- correlation_coefficient[R,R]
    #R_values <- R_values[lower.tri(R_values)]
    #R_values <- (sum(R_values, na.rm = TRUE) - sum(diag(R_values), na.rm = TRUE)) / 2
    #Rr2 <- sum(R_values,na.rm = TRUE)
    Rr2 <- sum(R_values, na.rm = TRUE) /2
    Rnna <- sum(!is.na(R_values))/2

    LR_values <- correlation_coefficient[L,R]
    LRr2 <- sum(LR_values,na.rm = TRUE)

    #OmegaD statistic
    #w = (choose(l,2) + choose(S-l,2))^(-1) * (Lr2 + Rr2) / (1/ l * (S-l)) * (LRr2)
    #w = (sum(!(is.na(Lr2))) + sum(!(is.na(Rr2))))^(-1) * (Lr2 + Rr2) / ((1/ sum(!(is.na(LRr2)))) * (LRr2))
    w = (Lnna + Rnna)^(-1) * (Lr2 + Rr2) / ((1/ sum(!(is.na(LR_values)))) * (LRr2))
  
    #Compiling each statistic into a table
    cat("Iteration: ",l,"\n")
    omega_values <- rbind(omega_values, c(l, w))
    #print(str(omega_values))
  
    rm(L_values, Lr2, LR_values, LRr2, R_values, Rr2)
    gc()
  }

  colnames(omega_values) <- c("l","omega")
  return(omega_values)
  
}
