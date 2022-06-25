# This script contains all code used before running fastStructure


library(tidyverse)


# Import the Data
Sys.setenv(VROOM_CONNECTION_SIZE="500000") # needed because the lines in this file are very long.
data.geno <- read_csv("../input/Rice_44K_genotypes.csv.gz",
                      na=c("NA","00")) #missing data is denoted as "NA" or "00"

data.geno <- data.geno %>% select(-`6_17160794...22253`) # remove a duplicate column
data.geno <- data.geno %>% rename(ID=...1) # Rename the first column to ID


# Construct PCA plot

# Convert data matrix to numbers
set.seed(0421)
data <- data.geno[,sample(2:ncol(data.geno), 10000)] # Subset the data
id <- select(data.geno, ID) # Extract the ID column
data.geno.10k <- data.frame(id, data) # Combine the ID column and 10k SNP subset
dim(data.geno.10k)