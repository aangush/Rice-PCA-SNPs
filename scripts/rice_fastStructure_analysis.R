# This script contains all code used to analyze the fastStructure output

# Below is the bash used to run fastStructure on Ubuntu 20.04 via WSL

# python structure.py -K 4 --input=/mnt/c/Users/aangu/bioinfo/bis180L/Rice-PCA-SNPs/output/rice.data.fastStructure.input 
# --output=/mnt/c/Users/aangu/bioinfo/bis180L/Rice-PCA-SNPs/output/rice.fastStructure.out --format=str


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


















