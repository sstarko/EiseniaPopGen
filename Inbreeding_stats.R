# Load required packages
library(adegenet)   # for gl2genind, handling genetic data
library(hierfstat)  # for basic.stats and genind2hierfstat conversion
library(boot)       # for bootstrapping
library(dplyr)      # for data manipulation
library(ggplot2)    # for plotting

gl_thinned
genind_obj <- gl2gi(gl_thinned)

# Make sure your genind object (genind_obj) has populations assigned.
# For example:
# pop_vector <- factor(c("MBS", "MBS", "MBD", "BA", "BAD", "SQS", "SQD", "ES", "ED", "T", "DS", "DD", ...))
# pop(genind_obj) <- pop_vector

# Get the unique populations from your genind object
pop_list <- levels(pop(genind_obj))

# Initialize vectors to store the mean stats and their CIs for each population.
n_pops <- length(pop_list)

popFis  <- numeric(n_pops)
popHo   <- numeric(n_pops)
popHe   <- numeric(n_pops)

lower_ci_Fis <- numeric(n_pops)
upper_ci_Fis <- numeric(n_pops)
lower_ci_Ho  <- numeric(n_pops)
upper_ci_Ho  <- numeric(n_pops)
lower_ci_He  <- numeric(n_pops)
upper_ci_He  <- numeric(n_pops)

# Function to compute the mean (used in bootstrapping)
boot_mean <- function(data, indices) {
  mean(data[indices], na.rm = TRUE)
}

# Loop over each population and compute the stats.
for(i in seq_along(pop_list)) {
  pop_name <- pop_list[i]
  
  # Subset the genind object for the current population
  subset_obj <- genind_obj[pop(genind_obj) == pop_name, ]
  
  # Convert to hierfstat format (population info becomes the first column)
  pop_data <- genind2hierfstat(subset_obj)
  
  # Compute basic statistics; basic.stats returns a list including Fis, Ho, and Hs
  bs_pop <- basic.stats(pop_data)
  
  # ---------------------------
  # Fis (Inbreeding coefficient)
  # ---------------------------
  fis_vec <- bs_pop$Fis  # vector of per-locus Fis for this population
  popFis[i] <- mean(fis_vec, na.rm = TRUE)
  fis_vec_clean <- fis_vec[!is.na(fis_vec)]
  
  if(length(fis_vec_clean) > 1) {
    boot_out <- boot(fis_vec_clean, statistic = boot_mean, R = 1000)
    ci_tmp <- boot.ci(boot_out, type = "perc")
    if(!is.null(ci_tmp$percent)) {
      lower_ci_Fis[i] <- ci_tmp$percent[4]
      upper_ci_Fis[i] <- ci_tmp$percent[5]
    } else {
      lower_ci_Fis[i] <- NA
      upper_ci_Fis[i] <- NA
    }
  } else {
    lower_ci_Fis[i] <- NA
    upper_ci_Fis[i] <- NA
  }
  
  # ---------------------------
  # Ho (Observed heterozygosity)
  # ---------------------------
  ho_vec <- bs_pop$Ho  # observed heterozygosity per locus
  popHo[i] <- mean(ho_vec, na.rm = TRUE)
  ho_vec_clean <- ho_vec[!is.na(ho_vec)]
  
  if(length(ho_vec_clean) > 1) {
    boot_out <- boot(ho_vec_clean, statistic = boot_mean, R = 1000)
    ci_tmp <- boot.ci(boot_out, type = "perc")
    if(!is.null(ci_tmp$percent)) {
      lower_ci_Ho[i] <- ci_tmp$percent[4]
      upper_ci_Ho[i] <- ci_tmp$percent[5]
    } else {
      lower_ci_Ho[i] <- NA
      upper_ci_Ho[i] <- NA
    }
  } else {
    lower_ci_Ho[i] <- NA
    upper_ci_Ho[i] <- NA
  }
  
  # ---------------------------
  # He (Expected heterozygosity)
  # ---------------------------
  # Note: hierfstat’s basic.stats returns expected heterozygosity as Hs.
  he_vec <- bs_pop$Hs  # expected heterozygosity per locus
  popHe[i] <- mean(he_vec, na.rm = TRUE)
  he_vec_clean <- he_vec[!is.na(he_vec)]
  
  if(length(he_vec_clean) > 1) {
    boot_out <- boot(he_vec_clean, statistic = boot_mean, R = 1000)
    ci_tmp <- boot.ci(boot_out, type = "perc")
    if(!is.null(ci_tmp$percent)) {
      lower_ci_He[i] <- ci_tmp$percent[4]
      upper_ci_He[i] <- ci_tmp$percent[5]
    } else {
      lower_ci_He[i] <- NA
      upper_ci_He[i] <- NA
    }
  } else {
    lower_ci_He[i] <- NA
    upper_ci_He[i] <- NA
  }
}

# Combine the results into one data frame
pop_stats_df <- data.frame(
  Population = pop_list,
  Fis = popFis,
  lower_ci_Fis = lower_ci_Fis,
  upper_ci_Fis = upper_ci_Fis,
  Ho = popHo,
  lower_ci_Ho = lower_ci_Ho,
  upper_ci_Ho = upper_ci_Ho,
  He = popHe,
  lower_ci_He = lower_ci_He,
  upper_ci_He = upper_ci_He,
  stringsAsFactors = FALSE
)

print(pop_stats_df)



# Define desired population order
pop_order <- c("MBS", "MBD", "BA", "BAD", "SQS", "SQD", "ES", "ED", "T", "DS", "DD")

# Ensure Population factor is ordered correctly
pop_stats_df <- pop_stats_df %>% 
  mutate(Population = factor(Population, levels = pop_order)) %>%
  arrange(Population)

# Plot mean Fis with error bars
plot_Fis<- ggplot(pop_stats_df, aes(x = Population, y = Fis)) +
  geom_col(fill = "grey") +
  geom_errorbar(aes(ymin = lower_ci_Fis, ymax = upper_ci_Fis), width = 0.2, color = "black") +
  scale_x_discrete(limits = pop_order) +
  theme_classic() +
  labs(
    title = "Inbreeding Coefficient (Fis)",
    x = "Population",
    y = "Inbreeding coefficient (Fis)"
  )+
  scale_y_continuous(limits = c(0, 0.8), expand = c(0, 0)) 



plot_Ho<- ggplot(pop_stats_df, aes(x = Population, y = Ho)) +
  geom_col(fill = "grey") +
  geom_errorbar(aes(ymin = lower_ci_Ho, ymax = upper_ci_Ho), width = 0.2, color = "black") +
  scale_x_discrete(limits = pop_order) +
  theme_classic() +
  labs(
    title = "Observed Heterozygosity (Ho)",
    x = "Population",
    y = "Mean observed heterozygosity (Ho)"
  )+
  scale_y_continuous(limits = c(0, 0.05),expand = c(0, 0)) 

plot_He <- ggplot(pop_stats_df, aes(x = Population, y = He)) +
  geom_col(fill = "grey") +
  geom_errorbar(aes(ymin = lower_ci_He, ymax = upper_ci_He), width = 0.2, color = "black") +
  scale_x_discrete(limits = pop_order) +
  theme_classic() +
  labs(
    title = "Expected Heterozygosity (He)",
    x = "Population",
    y = "Mean He (averaged across loci)"
  )+
  scale_y_continuous(limits = c(0, 0.06), expand = c(0, 0)) 

cowplot::plot_grid(plot_Fis, plot_Ho, plot_He, plot_Ne)


