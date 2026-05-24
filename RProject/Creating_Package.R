library("devtools")
library(roxygen2)

setwd("/home/eaa31/RProject/Correction_Package")
library(devtools);
load_all("."); # Working directory should be in the package SCC_R_package

Achaz_correction(dna_bin)
Omega_statistic(dna_bin)

gc()

library(roxygen2); # Read in the roxygen2 R package
roxygenise();      # Builds the help files
