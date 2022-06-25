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

# Fill missing data in with average genotype
geno.numeric.fill <-
  apply(geno.numeric, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm=T)
    x})

# Compute the PCs and return the first 10
geno.pca <- prcomp(geno.numeric.fill, 
                   rank.=10)

# Calculate the proportion of variance covered by each PC (for fun ;))
pcvar <- geno.pca$sdev^2 # square std dev to get variance
pcvar.pct <- tibble(pctvar=pcvar/sum(pcvar) * 100,
                    PC=1:length(pcvar))

# Plot the proportion of variance covered by each PC
pcvar.pct[1:10,] %>% 
  ggplot(mapping = aes(x=PC, y=pctvar)) +
  geom_col(fill="darkblue")

# Create a tibble with PCs and IDs
PCs <- as_tibble(geno.pca$x) %>% # The calculated PCs
  mutate(ID=data.geno.10k$ID) %>%
  select(ID, everything())

# Plot the PCs
PCs %>% 
  ggplot(mapping = aes(x=PC1, y=PC2)) +
  labs(title="PC1 vs PC2") +
  geom_point(color="darkblue")

PCs %>% 
  ggplot(mapping = aes(x=PC3, y=PC2)) +
  labs(title="PC3 vs PC2") +
  geom_point(color="darkblue")

# Create an MDS plot
# calculate the Euclidian distance between each rice variety
genDist <- as.matrix(dist(geno.numeric))

# perform multi-dimensional scaling
geno.mds <- as_tibble(cmdscale(genDist))

# Add the variety ID back into this
geno.mds$ID <- data.geno.10k$ID 

# Plot V1 vs V2 
geno.mds %>% 
  ggplot(mapping = aes(x=V1, y=V2)) +
  labs(title = "MDS Plot, V1 vs V2") +
  geom_point(color="darkblue") 

# Add phenotype data
data.pheno <- read_csv("../input/RiceDiversity.44K.MSU6.Phenotypes.csv")
data.pheno.pca <- inner_join(PCs, data.pheno, by=c("ID"="NSFTVID"))


# Color the PCA plots by different phenotypes
# Color PCA plot by Amylose content
data.pheno.pca %>%
  ggplot(mapping = aes(x=PC1, y=PC2, color=`Amylose content`)) +
  labs(title = "PC1 vs PC2 colored by Amylose content") +
  geom_point()

# Color PCA plot by pericarp color
data.pheno.pca %>% 
  ggplot(mapping = aes(x=PC1, y=PC2, color=`Pericarp color`)) +
  labs(title = "PC1 vs PC2 colored by Pericarp color") +  
  geom_point()

# Color PCA plot by Region
data.pheno.pca %>%
  ggplot(mapping = aes(x=PC1, y=PC2, color=`Region`)) +
  labs(title = "PC1 vs PC2 colored by Region") +
  scale_color_brewer(type="qual", palette = "Set1") +
  geom_point()



# fastStructure data preparation






















