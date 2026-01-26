library(vcfR)
library(StAMPP)

# Convert VCF to genotype matrix
genind_obj <- vcfR2genind(vcf_filtered2)  # Convert to genind
genind_obj@pop <- as.factor(sample_sites$site)
gl_fst <- dartR.base::gi2gl(genind_obj)


# Load required library
library(adegenet)

# Define the correct site order
ordered_sites <- c("MBS", "MBD", "BA", "BAD", "SQS", "SQD", "ES", "ED", "T", "DS", "DD")

# Extract population assignments from the genlight object
genlight_populations <- pop(gl_fst)  # This should be a factor

# Ensure the levels of populations match the desired order
genlight_populations <- factor(genlight_populations, levels = ordered_sites)

# Reorder the genlight object by these populations
gl_fst2 <- gl_fst[order(genlight_populations)]

# Convert to the correct format for StAMPP
stampp_data <- stamppConvert(gl_fst2, type = "genlight")

# Run pairwise Fst analysis on the ordered dataset
pairwise_fst <- stamppFst(stampp_data, nboots = 100, percent = 95)


# View results
print(pairwise_fst$Fsts)

library(ggplot2)
library(reshape2)

# Extract Fst matrix
fst_matrix <- as.matrix(pairwise_fst$Fsts)

# Convert to long format
fst_df <- melt(fst_matrix, na.rm = FALSE)
colnames(fst_df) <- c("Site1", "Site2", "Fst")

# Round values for display
fst_df$Fst_label <- sprintf("%.3f", fst_df$Fst)

# Keep only lower triangle (Site1 > Site2)
fst_df2 <- fst_df %>% drop_na()

# Define the correct site order
ordered_sites <- c("MBS", "MBD", "BA", "BAD", "SQS", "SQD", "ES", "ED", "T", "DS", "DD")

# Ensure only the lower triangle is used
#fst_df2 <- fst_df2[as.numeric(fst_df2$Site1) > as.numeric(fst_df2$Site2), ]

# Plot heatmap with fixed order
ggplot(fst_df2, aes(x = Site1, y = Site2, fill = Fst)) +
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(Fst), "", Fst_label)), color = "black", size = 3) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "#B30000",  
                       midpoint = quantile(fst_df2$Fst, 0.189, na.rm = TRUE),  
                       na.value = "white") +
  theme_classic() +
  labs(title = "", fill = "Fst") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(limits = ordered_sites) +  # Force correct order on X-axis
  scale_y_discrete(limits = c(ordered_sites))  # Reverse Y-axis for proper triangle
