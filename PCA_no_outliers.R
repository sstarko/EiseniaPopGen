snp_summary4 <- snp_summary3
snp_summary4$outlier <- rowSums(snp_summary3[, c("sigSST", "sigD", "sigSST_lfmm", "sigD_lfmm", "sigOUT")]) > 0

vcf_filtered2


###ADAPTIVE LOCI ONLY

# Get the SNP IDs that are flagged as outliers
outlier_ids <- snp_summary4$snpID[snp_summary4$outlier]

# Identify the row indices in the vcfR object that match these IDs.
# Adjust the column name ("ID") if your vcfR object stores locus names in a different column.
rows_to_keep <- which(vcf_filtered2@fix[,"ID"] %in% outlier_ids)

# Subset the vcfR object to retain only those rows.
vcf_filtered2_outliers <- vcf_filtered2[rows_to_keep, ]


#####PCA of just outlier SNPs

##Write VCF to file 
write.vcf(vcf_filtered2_outliers , "./intermediate_files/vcf_filtered2_outliers.vcf")

# Convert VCF to GDS format
vcf_file <- "./intermediate_files/vcf_filtered2_outliers.vcf"
gds_file <- "./intermediate_files/vcf_filtered2_outliers4.gds"
snpgdsVCF2GDS(vcf_file, gds_file, method="biallelic.only")

# Open the GDS file
genofile <- snpgdsOpen(gds_file)

# Run PCA
pca <- snpgdsPCA(genofile, autosome.only = FALSE)

# Get sample IDs and eigenvalues
sample.id <- pca$sample.id
pca_eigen <- pca$eigenval

# Create a data frame with the PCA results
pca_df <- data.frame(
  sample = sample.id,
  PC1 = pca$eigenvect[, 1],  # First principal component
  PC2 = pca$eigenvect[, 2]   # Second principal component
)

# Add the site information to the PCA data frame
pca_df$site <- extracted_sites[2:139]

pca_df$Depth <- pca_df$site %in% c("DD", "ED", "MBD", "SQD", "BAD") 

#FALSE INDICATES SHALLOW AND TRUE INDICATES DEEP


# Calculate percent variation explained by each PC
percent_var <- pca_eigen[1:32] / sum(pca_eigen[1:32]) * 100

# Format the axis labels with the percent variation
pc1_label <- paste0("PC1 (", round(percent_var[1], 2), "%)")
pc2_label <- paste0("PC2 (", round(percent_var[2], 2), "%)")

pca_df <- pca_df %>%
  mutate(site2 = case_when(
    site %in% c("ED", "ES") ~ "Ensenada",
    site %in% c("SQD", "SQS") ~ "SanQuentin",
    site %in% c("DD", "DS") ~ "Danvers",
    site == "T" ~ "Taylor",
    site %in% c("MBD", "MBS") ~ "Magdalena",
    site %in% c("BAD", "BA") ~ "BahiaAscension",
    TRUE ~ site # Retain original site name if not in the specified groups
  ))

library(geometry) # For convex hull calculation

# Calculate convex hulls
hull_data <- pca_df %>%
  group_by(site2) %>%
  slice(chull(PC1, PC2)) # Select points forming the convex hull



library(dplyr)
library(ggplot2)
library(geometry) # For convex hull calculation

# Calculate convex hulls
hull_data <- pca_df %>%
  group_by(site2) %>%
  slice(chull(PC1, PC2)) # Select points forming the convex hull


##HULLS
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  # Points colored by site
  geom_point(aes(color = site, shape = Depth), size = 2) +
  scale_color_manual(name = "Site Colors", values = c("red", "red", "blue", "blue", "lightpink", "lightpink",
                                                      "pink3", "pink3", "red4", "red4", "purple", "purple", "blue")) +
  # Convex hull polygons colored by site2
  geom_polygon(data = hull_data, aes(x = PC1, y = PC2, group = site2, fill = site2), 
               alpha = 0.2, color = "black") +
  scale_fill_manual(name = "Site2 Groups", values = c("darkred", "darkblue", "darkgreen", "gold", "orchid", "brown")) +
  # Theme and labels
  theme_classic() +
  labs(x = pc1_label, y = pc2_label)

##Ellipses
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  # Points colored by site
  geom_point(aes(color = site, shape = Depth), size = 2) +
  scale_color_manual(name = "Site Colors", values = c("red", "red", "blue4", "blue", "lightpink", "lightpink",
                                                      "pink3", "pink3", "red4", "red4", "purple", "purple", "blue"))+
  # Add a new color scale for ellipses
  new_scale_color() +
  # Ellipses colored by site2
  stat_ellipse(aes(color = site2, group = site2), level = 0.99) +
  scale_color_manual(name = "Site2 Ellipses", values = c("red", "blue", "lightpink", "pink3", "red4", "purple")) +
  # Theme and labels
  theme_classic() +
  labs(x = pc1_label, y = pc2_label)

##################
################NEUTRAL LOCI ONLY

# Get the SNP IDs that are flagged as outliers
outlier_ids <- snp_summary4$snpID[snp_summary4$outlier]

# Identify the row indices in the vcfR object that match these IDs.
# Adjust the column name ("ID") if your vcfR object stores locus names in a different column.
rows_to_keep <- which(!vcf_filtered2@fix[,"ID"] %in% outlier_ids)

# Subset the vcfR object to retain only those rows.
vcf_filtered2_outliers <- vcf_filtered2[rows_to_keep, ]


#####PCA of just outlier SNPs

##Write VCF to file 
write.vcf(vcf_filtered2_outliers , "./intermediate_files/vcf_filtered2_outliers.vcf")

# Convert VCF to GDS format
vcf_file <- "./intermediate_files/vcf_filtered2_outliers.vcf"
gds_file <- "./intermediate_files/vcf_filtered5_outliers.gds"
snpgdsVCF2GDS(vcf_file, gds_file, method="biallelic.only")

# Open the GDS file
genofile <- snpgdsOpen(gds_file)

# Run PCA
pca <- snpgdsPCA(genofile, autosome.only = FALSE)

# Get sample IDs and eigenvalues
sample.id <- pca$sample.id
pca_eigen <- pca$eigenval

# Create a data frame with the PCA results
pca_df <- data.frame(
  sample = sample.id,
  PC1 = pca$eigenvect[, 1],  # First principal component
  PC2 = pca$eigenvect[, 2]   # Second principal component
)

# Add the site information to the PCA data frame
pca_df$site <- extracted_sites[2:139]

pca_df$Depth <- pca_df$site %in% c("DD", "ED", "MBD", "SQD", "BAD") 

#FALSE INDICATES SHALLOW AND TRUE INDICATES DEEP


# Calculate percent variation explained by each PC
percent_var <- pca_eigen[1:32] / sum(pca_eigen[1:32]) * 100

# Format the axis labels with the percent variation
pc1_label <- paste0("PC1 (", round(percent_var[1], 2), "%)")
pc2_label <- paste0("PC2 (", round(percent_var[2], 2), "%)")

pca_df <- pca_df %>%
  mutate(site2 = case_when(
    site %in% c("ED", "ES") ~ "Ensenada",
    site %in% c("SQD", "SQS") ~ "SanQuentin",
    site %in% c("DD", "DS") ~ "Danvers",
    site == "T" ~ "Taylor",
    site %in% c("MBD", "MBS") ~ "Magdalena",
    site %in% c("BAD", "BA") ~ "BahiaAscension",
    TRUE ~ site # Retain original site name if not in the specified groups
  ))

library(geometry) # For convex hull calculation

# Calculate convex hulls
hull_data <- pca_df %>%
  group_by(site2) %>%
  slice(chull(PC1, PC2)) # Select points forming the convex hull



library(dplyr)
library(ggplot2)
library(geometry) # For convex hull calculation

# Calculate convex hulls
hull_data <- pca_df %>%
  group_by(site2) %>%
  slice(chull(PC1, PC2)) # Select points forming the convex hull


##HULLS
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  # Points colored by site
  geom_point(aes(color = site, shape = Depth), size = 2) +
  scale_color_manual(name = "Site Colors", values = c("red", "red", "blue", "blue", "lightpink", "lightpink",
                                                      "pink3", "pink3", "red4", "red4", "purple", "purple", "blue")) +
  # Convex hull polygons colored by site2
  geom_polygon(data = hull_data, aes(x = PC1, y = PC2, group = site2, fill = site2), 
               alpha = 0.2, color = "black") +
  scale_fill_manual(name = "Site2 Groups", values = c("darkred", "darkblue", "darkgreen", "gold", "orchid", "brown")) +
  # Theme and labels
  theme_classic() +
  labs(x = pc1_label, y = pc2_label)

##Ellipses
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  # Points colored by site
  geom_point(aes(color = site, shape = Depth), size = 2) +
  scale_color_manual(name = "Site Colors", values = c("red", "red", "blue4", "blue", "lightpink", "lightpink",
                                                      "pink3", "pink3", "red4", "red4", "purple", "purple", "blue"))+
  # Add a new color scale for ellipses
  new_scale_color() +
  # Ellipses colored by site2
  stat_ellipse(aes(color = site2, group = site2), level = 0.99) +
  scale_color_manual(name = "Site2 Ellipses", values = c("red", "blue", "lightpink", "pink3", "red4", "purple")) +
  # Theme and labels
  theme_classic() +
  labs(x = pc1_label, y = pc2_label)









###Pie charts showing some SNPs of interest
####TESTING LOCI OF INTERESTING

# Specify the locus of interest
locus_of_interest <- snp_summary3$snpID[640]

# Extract the specific locus from vcfR object

locus_data <- rownames(extract.gt(vcf_filtered2, element = "GT"))
to_keep <- which(locus_data %in% locus_of_interest)

vcf_interest<- vcf_filtered2[to_keep,]

# Extract gt

gt <- extract.gt(vcf_interest)

# Sample column names from vcf@gt
sample_names <- colnames(vcf_filtered2@gt)

# Use gsub with a regex pattern to extract the first 1-3 letters, stopping at the first number
extracted_sites <- gsub("^([A-Za-z]{1,3})\\d.*", "\\1", sample_names)

# Print the extracted site codes
print(extracted_sites)

# Create a function to calculate allele frequencies and plot pie charts
# Create a function to calculate allele frequencies and return data for plotting
plot_pie_for_sites <- function(genotype_data, site_names) {
  unique_sites <- unique(site_names)
  plot_data <- data.frame()
  
  for (site in unique_sites) {
    # Subset genotype data for the current site
    site_indices <- which(site_names == site)
    site_genotypes <- genotype_data[, site_indices]
    
    # Flatten the genotypes and count alleles
    allele_counts <- table(unlist(strsplit(site_genotypes, split = "[/|]")))
    total_alleles <- sum(allele_counts)
    allele_frequencies <- allele_counts / total_alleles
    
    # Create a data frame for plotting
    allele_df <- data.frame(
      Site = site,
      Allele = names(allele_frequencies),
      Frequency = as.numeric(allele_frequencies)
    )
    
    # Combine data for all sites
    plot_data <- bind_rows(plot_data, allele_df)
  }
  
  # Plot all pie charts as separate panels
  ggplot(plot_data, aes(x = "", y = Frequency, fill = Allele)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    theme_void() +
    scale_fill_manual(values = c("pink4", "lightblue"))+
    facet_wrap(~ Site) +
    theme(legend.position = "right") +
    labs(title = "Allele Frequencies by Site")
}

# Call the function with the extracted genotype data and site names
plot_pie_for_sites(gt, extracted_sites[-1])






