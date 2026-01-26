##Pie charts for figure

library(ggplot2)
library(dplyr)

# Define a layout mapping.
# Note: In Column 1 you have: MBS, BA, SQS, ES, T, DS
#       In Column 2 you have: MBD, BAD, SQD, ED, "" (blank), DS
layout_order <- data.frame(
  Site = c("MBS", "MBD",
           "BA",  "BAD",
           "SQS", "SQD",
           "ES",  "ED",
           "T",   "",     # blank for missing data
           "DS",  "DS"),
  Row = factor(rep(1:6, each = 2), levels = 6:1),  # levels reversed if you want row1 at the top
  Column = factor(rep(c("Col1", "Col2"), times = 6), levels = c("Col1", "Col2"))
)

# Revised function that accepts the layout mapping
plot_pie_for_sites <- function(genotype_data, site_names, layout_order) {
  # Compute allele frequencies only for sites with data
  unique_sites <- unique(site_names)
  plot_data <- data.frame()
  
  for (site in unique_sites) {
    # Find indices where this site occurs in your data
    site_indices <- which(site_names == site)
    site_genotypes <- genotype_data[, site_indices, drop = FALSE]
    
    # Count alleles from genotype strings (split by '/' or '|')
    allele_counts <- table(unlist(strsplit(site_genotypes, split = "[/|]")))
    total_alleles <- sum(allele_counts)
    allele_frequencies <- allele_counts / total_alleles
    
    # Create a temporary data frame for this site
    allele_df <- data.frame(
      Site = site,
      Allele = names(allele_frequencies),
      Frequency = as.numeric(allele_frequencies)
    )
    
    plot_data <- rbind(plot_data, allele_df)
  }
  
  # Merge the computed plot data with the layout mapping.
  # This ensures every row from layout_order is represented (even the blank one).
  plot_data_full <- layout_order %>%
    left_join(plot_data, by = "Site")
  
  # Create the pie charts using facet_grid to impose the 6 x 2 layout.
  ggplot(plot_data_full, aes(x = "", y = Frequency, fill = Allele)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    theme_void() +
    scale_fill_manual(values = c("pink4", "lightblue")) +
    facet_grid(Row ~ Column) +
    theme(legend.position = "right") +
    labs(title = "Allele Frequencies by Site")
}

# Example usage:

#8132
#7288
#4484
#1753

which(vcf_filtered2@fix[,3]=="104693:141")

vcf_interest <- vcf_filtered2[7319,]
gt <- extract.gt(vcf_interest)

# Assuming extracted_sites[-1] contains your site names,
# call the function with the layout_order:
plot_pie_for_sites(gt, extracted_sites[-1], layout_order)


#5844 depth
#7291 - possible depth?

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
