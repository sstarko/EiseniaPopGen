#This code calculates the error rates

#These are the eight duplicate pairs
#"BAD3_dup","BAD3" 
#"DS13_rep2","DS13_rep3" 
#"DS8_rep1","DS8_rep2" 
#"DS9_rep1","DS9_rep2" 
#"MBD3","MBD3_dup" 
#"MBS15" ,"MBS15_dup"
#"SQS14" ,"SQS14_dup" 
#"T5_dup","T5_rep2" 

#Create vcf with duplicates
vcf_filtered

# Load the vcfR package
library(vcfR)

# Define the duplicate sample names (ensure these match the column names in your vcf@gt)
dups <- c("BAD3_dup", "BAD3", 
          "DS13_rep2", "DS13_rep3", 
          "DS8_rep1", "DS8_rep2", 
          "DS9_rep1", "DS9_rep2", 
          "MBD3", "MBD3_dup", 
          "MBS15", "MBS15_dup",
          "SQS14", "SQS14_dup", 
          "T5_dup", "T5_rep2")

# Check available sample names (the first column in vcf@gt is "FORMAT")
available_samples <- colnames(vcf@gt)[-1]
if (!all(dups %in% available_samples)) {
  warning("Some of your duplicate sample names were not found in the VCF!")
}

# Subset the VCF object to include only the FORMAT column and your duplicate samples
vcf_dups <- vcf_filtered[, c("FORMAT", dups)]

# Optional: Verify that the new vcf object has only the desired samples
print(colnames(vcf_dups@gt))

# Assume 'vcf' is your vcfR object
# Extract genotype data (e.g., "GT" field)
gt <- extract.gt(vcf, element = "GT")

# Define a function to calculate error rate between two samples
calc_error_rate <- function(sample1, sample2) {
  # Identify positions where both samples have non-missing genotype calls
  valid <- !is.na(sample1) & !is.na(sample2)
  if(sum(valid) == 0) return(NA)  # Avoid division by zero if no valid comparisons
  
  # Count mismatches; note that for simple genotype strings (e.g., "0/0", "0/1", "1/1")
  # a direct string comparison is often sufficient.
  mismatches <- sum(sample1[valid] != sample2[valid])
  
  # Calculate error rate as the proportion of mismatches
  rate <- mismatches / sum(valid)
  return(rate)
}

# Example: Suppose your duplicate pairs are:
# Pair 1: "sampleA_rep1" and "sampleA_rep2"
# Pair 2: "sampleB_rep1" and "sampleB_rep2"
pairs <- list(c("BAD3_dup", "BAD3"),
              c("DS13_rep2", "DS13_rep3"),
              c("DS8_rep1", "DS8_rep2"), 
              c("DS9_rep1", "DS9_rep2"), 
              c("MBD3", "MBD3_dup"), 
              c("MBS15", "MBS15_dup"),
              c("SQS14", "SQS14_dup"), 
              c("T5_dup", "T5_rep2"))


 

# Calculate error rates for each duplicate pair
error_rates <- sapply(pairs, function(pair) {
  calc_error_rate(gt[, pair[1]], gt[, pair[2]])
})

# Display individual error rates
print(error_rates)

# Optionally, calculate the overall (mean) error rate across all duplicate pairs
overall_error_rate <- mean(error_rates, na.rm = TRUE)
print(overall_error_rate)

hist(error_rates, breaks = 9, las = 1, main = "Error rates")





# Assuming your VCF object is named 'vcf'
# Extract genotype calls (GT) and depth (DP) from the VCF
gt <- extract.gt(vcf, element = "GT")
dp <- extract.gt(vcf, element = "DP")

# Define duplicate sample names (in pairs)
dups <- c("BAD3_dup", "BAD3", 
          "DS13_rep2", "DS13_rep3", 
          "DS8_rep1", "DS8_rep2", 
          "DS9_rep1", "DS9_rep2", 
          "MBD3", "MBD3_dup", 
          "MBS15", "MBS15_dup",
          "SQS14", "SQS14_dup", 
          "T5_dup", "T5_rep2")

# Function to calculate missing data proportion and coverage summary for a pair of samples
compute_stats <- function(pair) {
  # Calculate missing data rates: missing if genotype is NA or "./."
  missing_rates <- sapply(pair, function(sample) {
    geno <- gt[, sample]
    mean(is.na(geno) | geno == "./.")
  })
  
  # Convert DP values to numeric and calculate summary statistics
  dp_numeric <- lapply(pair, function(sample) {
    as.numeric(dp[, sample])
  })
  coverage_stats <- lapply(dp_numeric, summary)
  
  list(missing = missing_rates, coverage = coverage_stats)
}

# Number of duplicate pairs (assuming dups length is even)
n_pairs <- length(dups) / 2

# Loop over each duplicate pair and print stats
for (i in 1:n_pairs) {
  pair <- dups[(2 * i - 1):(2 * i)]
  cat("\nDuplicate pair:", pair[1], "and", pair[2], "\n")
  
  stats <- compute_stats(pair)
  
  cat("Missing data proportion:\n")
  print(stats$missing)
  
  cat("Coverage summary:\n")
  print(stats$coverage)
}






# Assume 'vcf' is your vcfR object and 'dups' is defined as follows:
dups <- c("BAD3_dup", "BAD3", 
          "DS13_rep2", "DS13_rep3", 
          "DS8_rep1", "DS8_rep2", 
          "DS9_rep1", "DS9_rep2", 
          "MBD3", "MBD3_dup", 
          "MBS15", "MBS15_dup",
          "SQS14", "SQS14_dup", 
          "T5_dup", "T5_rep2")

# Extract genotype calls
gt <- extract.gt(vcf_dups, element = "GT")

# Function to calculate error rate between two samples
calc_error_rate <- function(sample1, sample2) {
  # Identify positions where both samples have valid genotype calls.
  valid <- !is.na(sample1) & !is.na(sample2) & sample1 != "./." & sample2 != "./."
  if (sum(valid) == 0) return(NA)  # Avoid division by zero
  mismatches <- sum(sample1[valid] != sample2[valid])
  mismatches / sum(valid)
}

# Initialize vectors to store error rates and missing data rates for each duplicate pair
n_pairs <- length(dups) / 2
error_rates <- numeric(n_pairs)
missing_rates <- numeric(n_pairs)

# Loop over each duplicate pair
for (i in 1:n_pairs) {
  pair <- dups[(2 * i - 1):(2 * i)]
  
  # Calculate error rate for the pair
  error_rates[i] <- calc_error_rate(gt[, pair[1]], gt[, pair[2]])
  
  # Calculate missing data rate for each sample in the pair (missing if NA or "./.")
  missing_rate1 <- mean(is.na(gt[, pair[1]]) | gt[, pair[1]] == "./.")
  missing_rate2 <- mean(is.na(gt[, pair[2]]) | gt[, pair[2]] == "./.")
  
  # Average the missing data rate for the pair
  missing_rates[i] <- mean(c(missing_rate1, missing_rate2))
}

# Create a data frame with the results
df <- data.frame(Pair = 1:n_pairs, ErrorRate = error_rates, MissingRate = missing_rates)
print(df)

# Plot Error Rate versus Missing Data Rate using ggplot2
ggplot(df, aes(x = MissingRate, y = ErrorRate)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, col = "blue") +
  labs(title = "Error Rate vs Missing Data Rate",
       x = "Average Missing Data Rate",
       y = "Error Rate")


