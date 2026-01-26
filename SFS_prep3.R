vcf_thinned


# Sample column names from vcf@gt
sample_names <- colnames(vcf_thinned@gt)

# Use gsub with a regex pattern to extract the first 1-3 letters, stopping at the first number
extracted_sites <- gsub("^([A-Za-z]{1,3})\\d.*", "\\1", sample_names)

# Extract genotype matrix
geno <- extract.gt(vcf_thinned, element = "GT")
colnames(geno) <- colnames(vcf_thinned@gt)[-1]  # Ensure sample names match

# Add site info
sample_sites <- data.frame(sample_name = colnames(geno), site = extracted_sites[-1])
rownames(sample_sites) <- sample_sites$sample_name
pop_conv <- data.frame(site = c("BA", "BAD", "MBD", "MBS", "SQS", "SQD", "ES", "ED", "DD", "DS", "T"),
                       cluster = c("1", "1", "1", "1", "2", "2", "3", "3", "4", "4", "4"),
                       group = c("S", "S", "S", "S", "S", "S", "S", "S", "N", "N", "N"))

sample_sites2 <- left_join(sample_sites, pop_conv)

library(vcfR)
library(dplyr)

# Load the metadata file (assuming it contains sample names and population info)
# Ensure this metadata dataframe has a column 'sample_name' and 'pop'
# meta <- read.csv("your_metadata_file.csv") 

# Count the number of individuals per population
pop_counts <- sample_sites2 %>% 
  group_by(site) %>%
  tally()

# Identify populations with at least 10 individuals
valid_pops <- pop_counts %>%
  filter(n >= 10) %>%
  pull(site)

# Get the samples that belong to those valid populations
valid_samples <- sample_sites2 %>%
  filter(site %in% valid_pops) %>%
  pull(sample_name)

# Match valid samples with those in the vcfR object
vcf_filtered_all <- vcf_thinned[, colnames(vcf_thinned@gt) %in% c("FORMAT", valid_samples)]



# Count the number of individuals per population
pop_counts <- sample_sites2 %>% 
  group_by(site) %>%
  tally()

# Identify populations with at least 10 individuals
valid_pops <- pop_counts %>%
  filter(n >= 10) %>%
  pull(site)

# Get the samples that belong to those valid populations
valid_samples1 <- sample_sites2 %>%
  filter(site == "BA") %>%
  pull(sample_name)

# Match valid samples with those in the vcfR object
vcf_filtered_BA <- vcf_thinned[, colnames(vcf_thinned@gt) %in% c("FORMAT", valid_samples1)]


# Get the samples that belong to those valid populations
valid_samples2 <- sample_sites2 %>%
  filter(site == "DS") %>%
  pull(sample_name)

# Match valid samples with those in the vcfR object
vcf_filtered_DS <- vcf_thinned[, colnames(vcf_thinned@gt) %in% c("FORMAT", valid_samples2)]


# Get the samples that belong to those valid populations
valid_samples3 <- sample_sites2 %>%
  filter(site == "ED") %>%
  pull(sample_name)

# Match valid samples with those in the vcfR object
vcf_filtered_ED <- vcf_thinned[, colnames(vcf_thinned@gt) %in% c("FORMAT", valid_samples3)]


# Get the samples that belong to those valid populations
valid_samples4 <- sample_sites2 %>%
  filter(site == "ES") %>%
  pull(sample_name)

# Match valid samples with those in the vcfR object
vcf_filtered_ES <- vcf_thinned[, colnames(vcf_thinned@gt) %in% c("FORMAT", valid_samples4)]

# Get the samples that belong to those valid populations
valid_samples5 <- sample_sites2 %>%
  filter(site == "MBS") %>%
  pull(sample_name)

# Match valid samples with those in the vcfR object
vcf_filtered_MBS <- vcf_thinned[, colnames(vcf_thinned@gt) %in% c("FORMAT", valid_samples5)]


# Get the samples that belong to those valid populations
valid_samples6 <- sample_sites2 %>%
  filter(site == "SQD") %>%
  pull(sample_name)

# Match valid samples with those in the vcfR object
vcf_filtered_SQD <- vcf_thinned[, colnames(vcf_thinned@gt) %in% c("FORMAT", valid_samples6)]

# Get the samples that belong to those valid populations
valid_samples7 <- sample_sites2 %>%
  filter(site == "SQS") %>%
  pull(sample_name)

# Match valid samples with those in the vcfR object
vcf_filtered_SQS <- vcf_thinned[, colnames(vcf_thinned@gt) %in% c("FORMAT", valid_samples7)]




# Ensure the directory exists
if (!dir.exists("./SFS_analysis2")) {
  dir.create("./SFS_analysis2", recursive = TRUE)
}

# Export the filtered VCF
vcfR::write.vcf(vcf_filtered_BA, file = "./SFS_analysis2/filtered_SFS_BA.vcf")
vcfR::write.vcf(vcf_filtered_DS, file = "./SFS_analysis2/filtered_SFS_DS.vcf")
vcfR::write.vcf(vcf_filtered_ED, file = "./SFS_analysis2/filtered_SFS_ED.vcf")
vcfR::write.vcf(vcf_filtered_ES, file = "./SFS_analysis2/filtered_SFS_ES.vcf")
vcfR::write.vcf(vcf_filtered_MBS, file = "./SFS_analysis2/filtered_SFS_MBS.vcf")
vcfR::write.vcf(vcf_filtered_SQD, file = "./SFS_analysis2/filtered_SFS_SQD.vcf")
vcfR::write.vcf(vcf_filtered_SQS, file = "./SFS_analysis2/filtered_SFS_SQS.vcf")
vcfR::write.vcf(vcf_filtered_all, file = "./SFS_analysis2/filtered_SFS_all.vcf")



# Define output file path
pops_file <- "./SFS_analysis2/pops_all.txt"


# Filter metadata to keep only valid samples
meta_filtered1 <- sample_sites2 %>%
  filter(site %in% valid_pops) %>%
  select(sample_name, site)  # Ensure the correct columns

# Write to pops.txt (tab-separated)
write.table(meta_filtered1, file = pops_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)




# Define output file path
pops_file <- "./SFS_analysis/pops_c1.txt"


# Filter metadata to keep only valid samples
meta_filtered1 <- sample_sites2 %>%
  filter(cluster == "1") %>%
  select(sample_name, cluster)  # Ensure the correct columns

# Write to pops.txt (tab-separated)
write.table(meta_filtered1, file = pops_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



# Define output file path
pops_file <- "./SFS_analysis/pops_c2.txt"


# Filter metadata to keep only valid samples
meta_filtered2 <- sample_sites2 %>%
  filter(cluster == "2") %>%
  select(sample_name, cluster)  # Ensure the correct columns

# Write to pops.txt (tab-separated)
write.table(meta_filtered2, file = pops_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# Define output file path
pops_file <- "./SFS_analysis/pops_c3.txt"


# Filter metadata to keep only valid samples
meta_filtered3 <- sample_sites2 %>%
  filter(cluster == "3") %>%
  select(sample_name, cluster)  # Ensure the correct columns

# Write to pops.txt (tab-separated)
write.table(meta_filtered3, file = pops_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# Define output file path
pops_file <- "./SFS_analysis/pops_c4.txt"


# Filter metadata to keep only valid samples
meta_filtered4 <- sample_sites2 %>%
  filter(cluster == "4") %>%
  select(sample_name, cluster)  # Ensure the correct columns

# Write to pops.txt (tab-separated)
write.table(meta_filtered4, file = pops_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


BA_samples <- sample_sites %>% filter(site == x)
BA_vcf <- vcf[, c("FORMAT", BA_samples$sample_name)]
vcfR::write.vcf(BA_vcf, file = "./BA_vcf.vcf")

library(pegas)

library(vcfR)
library(pegas)

# Convert to a genotype matrix
geno <- extract.gt(BA_vcf, element = "GT")

# Convert genotypes to allele counts (0,1,2 format)
geno_numeric <- apply(geno, 2, function(x) {
  as.numeric(gsub("[|/]", "", x))  # Replace "/" or "|" and convert to numeric
})

# Count derived allele frequencies per SNP
allele_counts <- apply(geno_numeric, 1, function(x) {
  sum(x, na.rm = TRUE)  # Count derived alleles
})

# Get the total number of alleles per SNP (2 per diploid individual)
n_alleles <- 2 * ncol(geno_numeric)

# Compute derived allele frequency per SNP
derived_freq <- allele_counts / n_alleles


# Convert allele frequencies into bins (for SFS)
sfs_bins <- table(round(derived_freq * length(derived_freq)))  # Bin frequencies

# Convert to a dataframe
sfs <- as.data.frame(sfs_bins)

# Save for use in Python
write.table(sfs, "sfs_input.txt", row.names = FALSE, col.names = FALSE)
