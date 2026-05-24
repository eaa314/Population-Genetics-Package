
Omega_statistic <- function(dna_bin) {

  data_test <- dna_bin
  colnames(data_test) <- paste0("snp",1:ncol(data_test))
  #data_test <- as.matrix(data_test)
  data_test <- as.character(data_test)

  #data_test <- read.dna("test_data.fasta",format = "fasta",
  #                      as.character = TRUE)

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

  correlation_coefficient <- cor(recoded_matrix, use = "pairwise.complete.obs")
  #correlation_coefficient[,1:10]

  summary(as.vector(correlation_coefficient[1:1000,1:1000]))
  #View(correlation_coefficient[31:50,31:50])
  #View(data_test[,"snp35"])

  correlation_coefficient <- correlation_coefficient^2
  diag(correlation_coefficient) <- 0
  snp_names <- colnames(recoded_matrix)
  S <- ncol(correlation_coefficient)
  omega_values <- data.frame(l = numeric(), w = numeric())

  gc()

  #for (l in 1:S) {
  for (l in 1:10) {
    #l <- 6
    L <- snp_names[1:l]
    R <- snp_names[-c(1:l)]

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
  
    cat("Iteration: ",l,"\n")
    omega_values <- rbind(omega_values, c(l, w))
    #print(str(omega_values))
  
    rm(L_values, Lr2, LR_values, LRr2, R_values, Rr2)
    gc()
  }

  print(omega_values)
  
}

Omega_statistic(dna_bin)
