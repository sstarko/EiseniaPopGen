vcf_filtered2

extracted_sites


# Load the libraries
library(vcfR)
library(adegenet)
library(zvau)

# Read your VCF file
vcf_filtered2

# Convert VCF to a genind object
genind_obj <- vcfR2genind(vcf_filtered2)

# Write out the genind object in Genepop format
writeGenPop(genind_obj, file.name = "./Genepop/all_samples.gen", comment = "All samples")


