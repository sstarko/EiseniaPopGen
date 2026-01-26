library(qvalue)
library(lfmm)

vcf_filtered2

# Sample column names from vcf@gt
sample_names <- colnames(vcf_filtered2@gt)

# Use gsub with a regex pattern to extract the first 1-3 letters, stopping at the first number
extracted_sites <- gsub("^([A-Za-z]{1,3})\\d.*", "\\1", sample_names)

# Print the extracted site codes
print(extracted_sites)

df <- data.frame(
  sites = extracted_sites,
  sample_names
)

df2 <- data.frame(
  sites = c("BA", "BAD", "SQS", "SQD", "ES", "ED", "MBS", "MBD", "DS", "DD", "T"),
  depth = c("shallow", "deep","shallow", "deep","shallow", "deep","shallow", "deep","shallow", "deep","shallow"),
  depth2 = c(0,1,0,1,0,1,0,1,0,1,0),
  depth3 = c(1, 10, 1, 10, 3, 10, 2, 12, 1, 1, 8)
)

meta<- full_join(df,df2)




#Determine K
# Convert to genind object for PCA
geno_data <- vcfR2genlight(vcf_filtered2)

# Run PCA
pca_result <- glPca(geno_data, nf = 10)  # Extracting first 10 principal components

# Plot the explained variance for each principal component
plot(pca_result$eig / sum(pca_result$eig), type = "b",
     xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     main = "PCA - Variance Explained")


# Extract genotypes and convert to numeric format (0, 1, 2)
convert_geno <- function(gt) {
  if (is.na(gt) || gt == "" || gt == "./.") {
    return(NA)
  } else if (gt == "0/0") {
    return(0)
  } else if (gt == "0/1" || gt == "1/0") {
    return(1)
  } else if (gt == "1/1") {
    return(2)
  } else {
    return(NA)
  }
}
geno_matrix <- apply(extract.gt(vcf_filtered2, element = "GT", as.numeric = FALSE), c(1, 2), convert_geno)

#Impute NAs as "0" which is the major allele and therefore the most common value (mode imputation)
geno_matrix[is.na(geno_matrix)] <- 0

# Save the genotype matrix in `.geno` format for LEA
write.geno(geno_matrix, "genotype_data.geno")
meta
meta2 <- meta[match(colnames(vcf_filtered2@gt)[2:139], meta$sample_names), ]


# Set the number of latent factors (K)
K <- 4 # Adjust based on your understanding of population structure


# Convert growth_rate to a dataframe
depth_df <- data.frame(depth = meta2$depth2[1:138])
depth2_df <- data.frame(depth = meta2$depth3[1:138])


## Fit an LFMM, i.e, compute B, U, V estimates
mod.lfmm <- lfmm_lasso(Y = t(geno_matrix), 
                       X = depth_df, 
                       K = 2,
                       nozero.prop = 0.01)

## performs association testing using the fitted model:
pv <- lfmm_test(Y = t(geno_matrix), 
                X = depth_df, 
                lfmm = mod.lfmm, 
                calibrate = "gif")


pvalues <- pv$pvalue 
p_values_adj <- p.adjust(pvalues, method = "fdr")
which(p_values_adj < 0.05)

## Check the p-values
hist(pv$pvalue, main="Unadjusted p-values")        

#hist(pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

# Let's change the GIF and readjust the p-values:
zscore <- pv$score[,1] 
gif <- pv$gif[1]     
new.gif1 <- 0.9 ### play around with this
adj.pv1 <- pchisq(zscore^2/new.gif1, df=1, lower = FALSE)
hist(pv$pvalue[,1], main="Unadjusted p-values")        
hist(pv$calibrated.pvalue[,1], main="GIF-adjusted p-values (GIF=1.12)")
hist(adj.pv1, main="REadjusted p-values (GIF=1)")

qv <- qvalue(adj.pv1)$qvalues
length(which(qv < 0.05)) ## h.w many SNPs have an FDR < 5% ?

qv2 <- qvalue(pv$calibrated.pvalue)$qvalues
length(which(qv2 < 0.1))

#This is the list of loci associated with depth according to LFMM
depth_loci_lfmm <- which(qv < 0.05) %>% names()
#depth_loci_lfmm <- which(qv < 0.05)

#imp_loci2 <- which(qv < 0.05) %>% names()
#causal_loci <- which(qv < 0.05) %>% names()


#####INTERACTION BETWEEN DEPTH AND TEMP

vcf_filtered2

# Sample column names from vcf@gt
sample_names <- colnames(vcf_filtered2@gt)

# Use gsub with a regex pattern to extract the first 1-3 letters, stopping at the first number
extracted_sites <- gsub("^([A-Za-z]{1,3})\\d.*", "\\1", sample_names)

# Print the extracted site codes
print(extracted_sites)

df <- data.frame(
  sites = extracted_sites,
  sample_names
)

df2 <- data.frame(
  sites = c("BA", "BAD", "SQS", "SQD", "ES", "ED", "MBS", "MBD", "DS", "DD", "T"),
  depth = c("shallow", "deep","shallow", "deep","shallow", "deep","shallow", "deep","shallow", "deep","shallow"),
  depth2 = c(0,1,0,1,0,1,0,1,0,1,0),
  depth3 = c(1, 10, 1, 10, 3, 10, 2, 12, 1, 1, 8)
)

meta<- full_join(df,df2)




#Determine K
# Convert to genind object for PCA
geno_data <- vcfR2genlight(vcf_filtered2)

# Run PCA
pca_result <- glPca(geno_data, nf = 10)  # Extracting first 10 principal components

# Plot the explained variance for each principal component
plot(pca_result$eig / sum(pca_result$eig), type = "b",
     xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     main = "PCA - Variance Explained")


# Extract genotypes and convert to numeric format (0, 1, 2)
convert_geno <- function(gt) {
  if (is.na(gt) || gt == "" || gt == "./.") {
    return(NA)
  } else if (gt == "0/0") {
    return(0)
  } else if (gt == "0/1" || gt == "1/0") {
    return(1)
  } else if (gt == "1/1") {
    return(2)
  } else {
    return(NA)
  }
}
geno_matrix <- apply(extract.gt(vcf_filtered2, element = "GT", as.numeric = FALSE), c(1, 2), convert_geno)

#Impute NAs as "0" which is the major allele and therefore the most common value (mode imputation)
geno_matrix[is.na(geno_matrix)] <- 0

# Save the genotype matrix in `.geno` format for LEA
write.geno(geno_matrix, "genotype_data.geno")
meta
meta2 <- meta[match(colnames(vcf_filtered2@gt)[2:139], meta$sample_names), ]


# Set the number of latent factors (K)
K <- 4 # Adjust based on your understanding of population structure


# Convert growth_rate to a dataframe
depth_df <- data.frame(depth = meta2$depth2[1:138])
depth2_df <- data.frame(depth = meta2$depth3[1:138])


## Fit an LFMM, i.e, compute B, U, V estimates
mod.lfmm <- lfmm_lasso(Y = t(geno_matrix), 
                       X = depth_df, 
                       K = 2,
                       nozero.prop = 0.01)

## performs association testing using the fitted model:
pv <- lfmm_test(Y = t(geno_matrix), 
                X = depth_df, 
                lfmm = mod.lfmm, 
                calibrate = "gif")


pvalues <- pv$pvalue 
p_values_adj <- p.adjust(pvalues, method = "fdr")
which(p_values_adj < 0.05)

## Check the p-values
hist(pv$pvalue, main="Unadjusted p-values")        

#hist(pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

# Let's change the GIF and readjust the p-values:
zscore <- pv$score[,1] 
gif <- pv$gif[1]     
new.gif1 <- 0.9 ### play around with this
adj.pv1 <- pchisq(zscore^2/new.gif1, df=1, lower = FALSE)
hist(pv$pvalue[,1], main="Unadjusted p-values")        
hist(pv$calibrated.pvalue[,1], main="GIF-adjusted p-values (GIF=1.12)")
hist(adj.pv1, main="REadjusted p-values (GIF=1)")

qv <- qvalue(adj.pv1)$qvalues
length(which(qv < 0.05)) ## h.w many SNPs have an FDR < 5% ?

qv2 <- qvalue(pv$calibrated.pvalue)$qvalues
length(which(qv2 < 0.1))

#This is the list of loci associated with depth according to LFMM
depth_loci_lfmm <- which(qv < 0.05) %>% names()
#depth_loci_lfmm <- which(qv < 0.05)

#imp_loci2 <- which(qv < 0.05) %>% names()
#causal_loci <- which(qv < 0.05) %>% names()







#install.packages("sdmpredictors")
#install.packages("raster")
library(sdmpredictors)
library(raster)
sst_layer <- load_layers("BO_sstmax")
sst_layer2 <- load_layers("BO_sstmean")

# Convert the dataframe to a SpatialPoints object
df3$lat <- as.numeric(df3$lat)
df3$long <- as.numeric(df3$long)
coordinates(df3) <- ~long + lat

# Extract SST data for the coordinates
df3$SST <- raster::extract(sst_layer, df3)
df3$SST2 <- raster::extract(sst_layer2, df3)

# Convert it back to a regular dataframe if needed
df3 <- as.data.frame(df3)

# View the updated dataframe
print(df3)


##Now look for outliers associated with temperature

## Fit an LFMM, i.e, compute B, U, V estimates
mod.lfmm <- lfmm_ridge(Y = t(geno_matrix), 
                       X = data.frame(df3$BO_sstmax), 
                       K = 4)

## performs association testing using the fitted model:
pv <- lfmm_test(Y = t(geno_matrix), 
                X = data.frame(df3$BO_sstmax), 
                lfmm = mod.lfmm, 
                calibrate = "gif")


pvalues <- pv$pvalue
p_values_adj <- p.adjust(pvalues, method = "fdr")
which(p_values_adj < 0.05)

## Check the p-values
hist(pv$pvalue, main="Unadjusted p-values")        

#hist(pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

# Let's change the GIF and readjust the p-values:
zscore <- pv$score[,1] 
gif <- pv$gif[1]     
new.gif1 <- 2 ### play around with this
adj.pv1 <- pchisq(zscore^2/new.gif1, df=1, lower = FALSE)
hist(pv$pvalue[,1], main="Unadjusted p-values")        
hist(pv$calibrated.pvalue[,1], main="GIF-adjusted p-values (GIF=1.12)")
hist(adj.pv1, main="REadjusted p-values (GIF=1)")

qv <- qvalue(adj.pv1)$qvalues
length(which(qv < 0.05)) ## h.w many SNPs have an FDR < 5% ?

qv2 <- qvalue(pv$calibrated.pvalue)$qvalues
length(which(qv2 < 0.1))

SST_loci_lfmm <- which(qv < 0.05) %>% names()
SST_loci_lfmm 
depth_loci_lfmm

snp_summary # Here is the data showing which SNPs are significant for RDA (both variables)

# Add new columns indicating whether each SNP is in the SST_loci_lfmm and depth_loci_lfmm lists
snp_summary2 <- snp_summary %>%
  mutate(
    sigSST_lfmm = snpID %in% SST_loci_lfmm,
    sigD_lfmm = snpID %in% depth_loci_lfmm
  )

# Check the result
print(snp_summary2, n = 10)  # Print the first 10 rows


# Load necessary library
library(VennDiagram)


# Load necessary library
library(VennDiagram)
library(grid)

# Example: Load your dataset (Replace with your actual file path)
df <- snp_summary2

# Ensure logical values (Convert from character if needed)
df$sigSST <- as.logical(df$sigSST)
df$sigD <- as.logical(df$sigD)
df$sigSST_lfmm <- as.logical(df$sigSST_lfmm)
df$sigD_lfmm <- as.logical(df$sigD_lfmm)

# Identify SNPs that are TRUE for each category
sigSST_set <- df$snpID[df$sigSST == TRUE]
sigD_set <- df$snpID[df$sigD == TRUE]
sigSST_lfmm_set <- df$snpID[df$sigSST_lfmm == TRUE]
sigD_lfmm_set <- df$snpID[df$sigD_lfmm == TRUE]

# Create a 4-set Venn diagram
venn.plot <- draw.quad.venn(
  area1 = length(sigSST_set),
  area2 = length(sigD_set),
  area3 = length(sigSST_lfmm_set),
  area4 = length(sigD_lfmm_set),
  n12 = length(intersect(sigSST_set, sigD_set)),
  n13 = length(intersect(sigSST_set, sigSST_lfmm_set)),
  n14 = length(intersect(sigSST_set, sigD_lfmm_set)),
  n23 = length(intersect(sigD_set, sigSST_lfmm_set)),
  n24 = length(intersect(sigD_set, sigD_lfmm_set)),
  n34 = length(intersect(sigSST_lfmm_set, sigD_lfmm_set)),
  n123 = length(intersect(intersect(sigSST_set, sigD_set), sigSST_lfmm_set)),
  n124 = length(intersect(intersect(sigSST_set, sigD_set), sigD_lfmm_set)),
  n134 = length(intersect(intersect(sigSST_set, sigSST_lfmm_set), sigD_lfmm_set)),
  n234 = length(intersect(intersect(sigD_set, sigSST_lfmm_set), sigD_lfmm_set)),
  n1234 = length(intersect(intersect(intersect(sigSST_set, sigD_set), sigSST_lfmm_set), sigD_lfmm_set)),
  category = c("sigSST", "sigD", "sigSST_lfmm", "sigD_lfmm"),
  fill = c("red", "blue", "green", "purple"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.2
)


# Display the plot
plot.new()
grid.draw(venn.plot)





####TESTING LOCI OF INTERESTING
causal_loci

# Specify the locus of interest
locus_of_interest <- causal_loci[14]

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
    facet_wrap(~ Site) +
    theme(legend.position = "right") +
    labs(title = "Allele Frequencies by Site")
}

# Call the function with the extracted genotype data and site names
plot_pie_for_sites(gt, extracted_sites[-1])


##Coords

coordinates <- data.frame(
  site = c("BA", "BAD", "DD", "DS", "ED", "ES", "MBD", "MBS", "SQD", "SQS", "T"),
  lat = c(
    "27.139796",
    "27.139796",
    "48.876851", 
    "48.876851", 
    "31.867667", 
    "31.867667", 
    "24.647856", 
    "24.647856",
    "30.567708", 
    "30.567708", 
    "48.827562"
  ),
  long = c(
    "-114.291897",
    "-114.291897",
    "-125.092248",
    "-125.092248",
    "-116.680746",
    "-116.680746",
    "-111.880664",
    "-111.880664",
    "-116.036209",
    "-116.036209",
    "-125.197424"
  )
)


df <- data.frame(sample = sample_names,
site = extracted_sites)

df2 <- left_join(df, coordinates, by= "site")
df3 <- df2[2:139,]


#install.packages("sdmpredictors")
#install.packages("raster")
library(sdmpredictors)
library(raster)
sst_layer <- load_layers("BO_sstmax")
sst_layer2 <- load_layers("BO_sstmean")

# Convert the dataframe to a SpatialPoints object
df3$lat <- as.numeric(df3$lat)
df3$long <- as.numeric(df3$long)
coordinates(df3) <- ~long + lat

# Extract SST data for the coordinates
df3$SST <- raster::extract(sst_layer, df3)
df3$SST2 <- raster::extract(sst_layer2, df3)

# Convert it back to a regular dataframe if needed
df3 <- as.data.frame(df3)

# View the updated dataframe
print(df3)


##Now look for outliers associated with temperature

## Fit an LFMM, i.e, compute B, U, V estimates
mod.lfmm <- lfmm_ridge(Y = t(geno_matrix), 
                       X = data.frame(df3$BO_sstmax), 
                       K = 4)

## performs association testing using the fitted model:
pv <- lfmm_test(Y = t(geno_matrix), 
                X = data.frame(df3$BO_sstmax), 
                lfmm = mod.lfmm, 
                calibrate = "gif")


pvalues <- pv$calibrated.pvalue 
p_values_adj <- p.adjust(pvalues, method = "fdr")
which(p_values_adj < 0.01)

## Check the p-values
hist(pv$pvalue, main="Unadjusted p-values")        

hist(pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")


vcf_interest<- vcf_filtered2[5401,]

gt <- extract.gt(vcf_interest)
# Call the function with the extracted genotype data and site names
plot_pie_for_sites(gt, extracted_sites[-1])


vcf_filtered2

#8038 is a cool one that seems to be associated with both temp and depth
#7618 opposite


#For deep
#944?
#4478?
#5401?






#####PRINT PIE CHARTS OF EACH OUTLIER

# Create the subdirectory "pie_charts" if it doesn't already exist
if (!dir.exists("pie_charts")) {
  dir.create("pie_charts")
}

# Loop over each SNP (assuming each row in vcf_filtered2 is one SNP)
for(i in 1:nrow(vcf_filtered2_outliers)) {
  
  # Extract the SNP of interest
  vcf_interest <- vcf_filtered2_outliers[i, ]
  gt <- extract.gt(vcf_interest)
  
  # Define the file name for the plot
  file_name <- paste0("pie_charts/pie_chart_SNP_", i, ".png")
  
  # Open a PNG device (adjust width/height as needed)
  png(filename = file_name, width = 800, height = 600)
  
  # If plot_pie_for_sites returns a ggplot object, assign it to p and explicitly print it
  p <- plot_pie_for_sites(gt, extracted_sites[-1])
  print(p)
  
  # Close the device, saving the plot
  dev.off()
}
