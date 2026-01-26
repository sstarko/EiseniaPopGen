##Thin the dataset to include first SNP per locus
library(dartR.base)
gl1 <- vcfR2genlight(vcf_filtered2)

# Extract the RADtag from the marker names
radtag <- sapply(strsplit(locNames(gl1), ":"), `[`, 1)

# Identify the first SNP encountered for each RADtag
keep <- !duplicated(radtag)

# Subset gl1 to keep only one SNP per RADtag
gl1_thinned <- gl1[, keep]

gl2vcf(gl1_thinned)


fix_df <- data.frame(vcf_filtered2@fix)

# Group by RADtag and select the first SNP per RADtag.
fix_df_thinned <- fix_df %>% group_by(CHROM) %>% slice(1)


# Get the row numbers corresponding to the selected SNPs.
selected_rows <- which(fix_df$ID %in% fix_df_thinned$ID)

# Subset the original VCF object.
vcf_thinned <- vcf_filtered2[selected_rows, ]

gl_thinned <- vcfR2genlight(vcf_thinned)
