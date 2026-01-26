###RDA analysis (and outlier detection)
#Load libraries
library(adegenet)
library(dartR.base)
library(vegan)   # For RDA
library(adegenet) # For working with genind objects

#Check that data are loaded from a previous script
gl_obj2
#Rename
gl.final <- gl_obj2

#Impute missing data

genind_obj <- gl2gi(gl.final)
genind_obj <- scaleGen(genind_obj, NA.method = "mean")

#Check that metadata are loaded
df3

# Extract SNP matrix from genind object
snp_data <- tab(as.genind(genind_obj), NA.method = "mean")




# Fit RDA
rda_model <- rda(snp_data ~ df3$BO_sstmax + depth_df$depth)
summary(rda_model)

# Extract SNP loadings for the first RDA axis
loadings <- scores(rda_model, choices = 1, display = "species")

# Identify outliers as those SNPs with loadings > 2 standard deviations from the mean
outlier_threshold <- 2.5 * sd(loadings)
outliers <- which(abs(loadings - mean(loadings)) > outlier_threshold)

# View the outlier SNPs
outlier_snps <- colnames(snp_data)[outliers]
outlier_snps

# Extract site scores (samples) and SNP scores (species)
site_scores <- as.data.frame(scores(rda_model, display = "sites"))
site_scores$site <- df3$site  # Add site labels from meta

snp_scores <- as.data.frame(scores(rda_model, display = "species"))
snp_scores$SNP <- colnames(snp_data)  # Add SNP names

# Identify significant SNPs on the first RDA axis using 2 SDs from the mean
snp_scores$significant <- abs(snp_scores$RDA1 - mean(snp_scores$RDA1)) > 2.5 * sd(snp_scores$RDA1)

# Identify significant SNPs on the second RDA axis using 2 SDs from the mean
snp_scores$significant2 <- abs(snp_scores$RDA2 - mean(snp_scores$RDA2)) > 2.5 * sd(snp_scores$RDA2)

# Calculate z-scores for RDA1 loadings
z_scores <- (snp_scores$RDA1 - mean(snp_scores$RDA1)) / sd(snp_scores$RDA1)
snp_scores$p_value_RDA1 <- 2 * (1 - pnorm(abs(z_scores)))  # Two-tailed p-value

# Similarly for RDA2 loadings
z_scores2 <- (snp_scores$RDA2 - mean(snp_scores$RDA2)) / sd(snp_scores$RDA2)
snp_scores$p_value_RDA2 <- 2 * (1 - pnorm(abs(z_scores2)))

qvalue(snp_scores$p_value_RDA1)$qvalue %>% min()



###################################
############THRESHOLD METHOD#######

#Fit RDA
rda_model <- rda(snp_data ~ df3$BO_sstmax + depth_df$depth)
summary(rda_model2)


# Extract SNP loadings for the first RDA axis
loadings <- scores(rda_model, choices = 1, display = "species")

# Identify outliers in the top 5% of loadings
threshold <- quantile(abs(loadings), 0.95)
outliers <- which(abs(loadings) > threshold)

# View the outlier SNPs
outlier_snps <- colnames(snp_data)[outliers]
outlier_snps


# Extract site scores (samples) and SNP scores (species)
site_scores <- as.data.frame(scores(rda_model, display = "sites"))
site_scores$site <- df3$site  # Add site labels from meta

snp_scores <- as.data.frame(scores(rda_model, display = "species"))
snp_scores$SNP <- colnames(snp_data)  # Add SNP names

# Identify significant SNPs based on the top 5% loadings on the first RDA axis
loading_threshold <- quantile(abs(snp_scores$RDA1), 0.99)
snp_scores$significant <- abs(snp_scores$RDA1) > loading_threshold  # Logical vector for significant SNPs

# Identify significant SNPs based on the top 5% loadings on the first RDA axis
loading_threshold2 <- quantile(abs(snp_scores$RDA2), 0.99)
snp_scores$significant2 <- abs(snp_scores$RDA2) > loading_threshold2  # Logical vector for significant SNPs

#############################

# Plotting samples and SNPs together on the RDA
ggplot() +
  # Plot samples (sites) with color by 'site'
  geom_point(data = site_scores, aes(x = RDA1, y = RDA2, color = site, )) +
  #geom_text(data = site_scores, aes(x = RDA1, y = RDA2, label = site), vjust = -1, size = 3) +
  
  # Plot SNPs, highlighting significant ones
  geom_point(data = snp_scores, aes(x = RDA1, y = RDA2), color = "grey", alpha = 0.5) +
  geom_point(data = filter(snp_scores, significant == TRUE), aes(x = RDA1, y = RDA2), color = "red") +
  geom_point(data = filter(snp_scores, significant2 == TRUE), aes(x = RDA1, y = RDA2), color = "pink") +
  
  # Labels and themes
  labs(x = "RDA1 (11.1 %)", y = "RDA2 (1.3%)") +
  theme_classic() +
  theme(legend.position = "bottom")

#Plotting just the SNPs
ggplot() +
  # Plot samples (sites) with color by 'site'
  #geom_text(data = site_scores, aes(x = RDA1, y = RDA2, label = site), vjust = -1, size = 3) +
  
  # Plot SNPs, highlighting significant ones
  geom_point(data = snp_scores, aes(x = RDA1, y = RDA2), color = "grey", alpha = 0.5) +
  geom_point(data = filter(snp_scores, significant == TRUE), aes(x = RDA1, y = RDA2), color = "red") +
  geom_point(data = filter(snp_scores, significant2 == TRUE), aes(x = RDA1, y = RDA2), color = "pink") +
  
  # Labels and themes
  labs(x = "RDA1 (11.1 %)", y = "RDA2 (1.3%)") +
  theme_classic() +
  theme(legend.position = "bottom")

which(gl.final@loc.names == "21579:23")
which(gl.final@loc.names == "86006:101")
which(gl.final@loc.names == "36293:18")


#(All from LFMM)
#1753 is a cool snp
#8038 has deep clustering with south
#962 also has deep clustering with south
#7618 also has deep clustering with south
#3056 also has deep clustering with south
#2909 looks like it is variable in mid-range only, more abundant at deep sites
#282 Appears under selection in Ensenada
#4119 deep appears to cluster with south (at San Quintin)
#4117 appears under selection in San Quintin

#6060 shallow clusters with south
#8132
#401 is an interesting one for under selection across

vcf_interest<- vcf_filtered2[3056,]

gt <- extract.gt(vcf_interest)
# Call the function with the extracted genotype data and site names
plot_pie_for_sites(gt, extracted_sites[-1])



snp_scores2 <- snp_scores
snp_scores2$snpID <- rownames(snp_scores2)
snp_scores2$snpID <- gsub(".T", "", snp_scores2$snpID)
snp_scores2$snpID <- gsub(".G", "", snp_scores2$snpID)
snp_scores2$snpID <- gsub(".A", "", snp_scores2$snpID)
snp_scores2$snpID <- gsub(".C", "", snp_scores2$snpID)

snp_summary <- snp_scores2 %>% group_by(snpID) %>% 
  reframe(sigSST = unique(significant) , sigD = unique(significant2))





rda_model_all <- rda(snp_data ~ df3$BO_sstmax + depth_df$depth)
site_scores_all <- as.data.frame(scores(rda_model_all, display = "sites"))
site_scores_all$site <- df3$site  # Add site labels




# Filter SNP data to only include significant SNPs
sig_snps <- colnames(snp_data)[snp_scores$significant | snp_scores$significant2]  # SNPs with high RDA1 or RDA2 loadings
snp_data_sig <- snp_data[, sig_snps, drop = FALSE]  # Keep only significant SNPs

# Run RDA on the reduced SNP dataset
rda_model_sig <- rda(snp_data_sig ~ df3$BO_sstmax + depth_df$depth)
site_scores_sig <- as.data.frame(scores(rda_model_sig, display = "sites"))
site_scores_sig$site <- df3$site  # Add site labels


# Add shape column based on whether site name contains "D"
site_scores_all$shape <- ifelse(grepl("D", site_scores_all$site), 15, 16)  # 15 = square, 16 = circle
site_scores_sig$shape <- ifelse(grepl("D", site_scores_sig$site), 15, 16)


# Extract proportion of variance explained by each axis
rda_var_all <- summary(rda_model_all)$cont$importance[2, 1:2] * 100  # Convert to percentage
rda_var_sig <- summary(rda_model_sig)$cont$importance[2, 1:2] * 100


ggplot(site_scores_all, aes(x = RDA1, y = RDA2, color = site, shape = as.factor(shape))) +
  geom_point(size = 3) +
  scale_shape_manual(values = c("15" = 16, "16" = 17)) +  # Define square and circle shapes
  labs(
    title = "RDA - All SNPs",
    x = paste0("RDA1 (", round(rda_var_all[1], 1), "%)"),
    y = paste0("RDA2 (", round(rda_var_all[2], 1), "%)")
  ) +
  theme_classic() +
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("pink3", "pink3", "blue4", "blue", "lightpink", "lightpink","red", "red",
                                 "red4", "red4", "purple", "purple","blue"))

ggplot(site_scores_sig, aes(x = RDA1, y = RDA2, color = site, shape = as.factor(shape))) +
  geom_point(size = 3) +
  scale_shape_manual(values = c("15" = 15, "16" = 16)) +  # Define square and circle shapes
  labs(
    title = "RDA - Significant SNPs Only",
    x = paste0("RDA1 (", round(rda_var_sig[1], 1), "%)"),
    y = paste0("RDA2 (", round(rda_var_sig[2], 1), "%)")
  ) +
  theme_classic() +
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("pink3", "pink3", "blue4", "blue", "lightpink", "lightpink","red", "red",
                                "red4", "red4", "purple", "purple","blue"))



snp_summary # Here is the data showing which SNPs are significant for RDA (both variables)
