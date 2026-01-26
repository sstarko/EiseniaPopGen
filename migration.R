library(adegenet)
library(diveRsity)
library(DNAprofiles)
library(divMigrateR)

library(adegenet)
library(SNPRelate)
library(gdsfmt)

#genlight_obj <- gl_obj2

#Remove BAD, DD, MBD due to low sample size (under 10)
# Define populations to remove
remove_pops <- c("BAD", "DD", "MBD", "T")

# Get population info
pop_info <- pop(gl_obj2)

# Identify individuals that belong to populations NOT in the remove list
keep_inds <- !(pop_info %in% remove_pops)

# Filter the genlight object
genl_filtered <- gl_obj2[keep_inds, ]

# Check population levels after filtering
pop(genl_filtered) <- droplevels(pop(genl_filtered))

# Print summary to confirm
print(genl_filtered)

genl_filtered@pop %>% table()
genlight_obj <- genl_filtered

library(adegenet)

# Ensure population information exists
if (is.null(pop(genlight_obj))) {
  stop("Error: Population data is missing in your genlight object! Assign populations before converting.")
}

# Extract genotype matrix (coded as 0,1,2)
geno_matrix <- as.matrix(genlight_obj)

# Extract sample names
sample_names <- indNames(genlight_obj)

# Extract SNP locus names
loci_names <- locNames(genlight_obj)

# Extract population information
populations <- as.character(pop(genlight_obj))

# Define output file
genepop_filename <- "eisenia_migration.gen"

# Open file for writing
fileConn <- file(genepop_filename, "w")

# Write Genepop header
writeLines("Genepop file for divMigrate", fileConn)

# Write SNP locus names
writeLines(paste(loci_names, collapse = ","), fileConn)

# Group individuals by population
unique_pops <- unique(populations)
for (pop_id in unique_pops) {
  writeLines("POP", fileConn)  # Genepop format requires "POP" before each population
  
  pop_indices <- which(populations == pop_id)
  
  for (i in pop_indices) {
    # Convert SNP genotypes to Genepop format (two-digit alleles)
    geno_formatted <- apply(geno_matrix[i, , drop = FALSE], 2, function(x) {
      if (is.na(x)) return("0000")  # Missing data
      if (x == 0) return("0101")  # Homozygous for first allele
      if (x == 1) return("0102")  # Heterozygous
      if (x == 2) return("0202")  # Homozygous for second allele
      return("0000")  # Default in case of issues
    })
    
    # Write individual data
    writeLines(paste(sample_names[i], ",", paste(geno_formatted, collapse = " ")), fileConn)
  }
}

# Close the file
close(fileConn)

library(diveRsity)

# Run divMigrate with the properly formatted Genepop file
migrate_results <- divMigrate(infile = "eisenia_migration.gen", plot_network = TRUE, boots = 1000, filter_threshold = 0.05)
migrate_results <- divMigrate(infile = "eisenia_migration.gen", plot_network = TRUE, filter_threshold = 0.05)




