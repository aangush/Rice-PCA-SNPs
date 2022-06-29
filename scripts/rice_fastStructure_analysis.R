# This script contains all code used to analyze the fastStructure output

# Below is the bash used to run fastStructure on Ubuntu 20.04 via WSL

# python structure.py -K 4 --input=/mnt/c/Users/aangu/bioinfo/bis180L/Rice-PCA-SNPs/output/rice.data.fastStructure.input 
# --output=/mnt/c/Users/aangu/bioinfo/bis180L/Rice-PCA-SNPs/output/rice.fastStructure.out --format=str


library(tidyverse)

# Load the fastStructure results
fs_results <- read_delim("../output/rice.fastStructure.out.4.meanQ", delim="  ", col_names = FALSE, col_types = 'nnnn')

# Add sample IDs back and assign population names
fs_results <- fs_results %>% 
  mutate(ID=data.geno.10k$ID) %>% 
  select(ID, pop1=X1, pop2=X2, pop3=X3, pop4=X4)

# Assign each individual to a sub pop.
fs_results$assignedPop <- apply(fs_results[,-1], 1, which.max)


# Order samples based on population identity
fs_results$maxPr <- apply(fs_results[,2:5],1,max) 
fs_results <- fs_results %>% 
  arrange(assignedPop,desc(maxPr)) %>%
  mutate(plot.order=row_number())

# Convert data to long format
fs_results_long <- fs_results %>% pivot_longer(pop1:pop4, 
                                               names_to="population",
                                               values_to="proportion")
# Plot the fastStructure results
fs_results_long %>%
  ggplot(aes(x=plot.order, y=proportion, color=population, fill=population)) + 
  geom_col()  +
  ylab("genome proportion") + 
  scale_color_brewer(type="div") + scale_fill_brewer(type="div")


# Compare sub pop. assignments to PCA plot

# Combine PCs with fastStructure data
fs_results <- fs_results %>% mutate(assignedPop=as.character(assignedPop))
geno.pca.pop <- inner_join(fs_results, PCs, by="ID")


# PCA plots colored by pop. assignment
# PC1 vs PC2
geno.pca.pop %>% ggplot(mapping = aes(x=PC1, y=PC2, color=assignedPop)) +
  labs(title="PC1 vs PC2 colored by sub pop. assignment") +
  geom_point()

# PC2 vs PC3
geno.pca.pop %>% ggplot(mapping = aes(x=PC3, y=PC2, color=assignedPop)) +
  labs(title="PC3 vs PC2 colored by sub pop. assignment") +
  geom_point()


