library(snmf)
library(LEA)
library(dartR.base)

##SNMF
############################################
##################ALL SNPS##################
############################################

#Genlight from previous file

gl_obj2  #138 samples, 9,315 variants


meta <- data.frame(site = c("BA", "BAD", "DD", "DS", "ED", "ES", "MBD", "MBS", "SQD", "SQS", "T"), 
                   site2 = c("BA", "BA","BS", "BS", "E", "E", "MB", "MB", "SQ", "SQ", "BS"),
                   depth = c("Shallow", "Deep", "Deep", "Shallow", "Deep", "Shallow", "Deep", "Shallow", "Deep", "Shallow", "Shallow")
)


meta


# Sample column names from vcf@gt
sample_names <- colnames(vcf_filtered2@gt)

# Use gsub with a regex pattern to extract the first 1-3 letters, stopping at the first number
site <- gsub("^([A-Za-z]{1,3})\\d.*", "\\1", sample_names)[2:139]

meta2 <- data.frame(sample = sample_names[2:139], site = site)

meta.data <- left_join(meta, meta2, by = "site")

meta.data

gl_obj2

gl2faststructure(gl_obj2, outfile = "./gl_structure.fstr",
                 outpath = './intermediate_files2/')


LEA::struct2geno("./intermediate_files2/gl_structure.fstr", ploidy = 2, FORMAT = 2)
obj.snmf = snmf("./intermediate_files2/gl_structure.fstr.geno", 1:10,
                entropy = TRUE, ploidy = 2,project = "new",  repetitions = 10, alpha = 100)
qmatrix = Q(obj.snmf, K = 5, run = 10)

barplot(t(qmatrix), col = c("#D1E5F0","#D6604D","#FDDBC7","#F4A582","#92C5DE","#4393C3"), border = NA, space = 0,
        xlab = "Individuals", ylab = "Admixture coefficients")

plot(obj.snmf, pch = 19, cex = 1.5)


# Ensure meta.data has the correct mapping of samples and sites
meta.data <- left_join(meta2, meta, by = "site")

# Reorder meta.data by the desired site order
desired_order <- c("DD", "DS","T", "ED", "ES",  "SQD", "SQS", "BA", "BAD", "MBD", "MBS")
meta.data <- meta.data %>%
  mutate(site = factor(site, levels = desired_order)) %>%
  arrange(site)


# Reorder qmatrix rows to match the order of samples in meta.data
qmatrix2 <- qmatrix[match(meta.data$sample, indNames(gl_obj2)), ]

# Create the barplot with the reordered qmatrix
all_snp_barplot_k4<- barplot(t(qmatrix2), col = c("#4393C3","#D6604D","#D1E5F0",   "#FDDBC7"),
        border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients")+
  axis(1, at = seq_along(meta.data$site), labels = meta.data$site, las = 2, cex.axis = 0.7)

#Add site labels below the barplot
axis(1, at = seq_along(meta.data$site), labels = meta.data$site, las = 2, cex.axis = 0.7)


# Create the barplot with the reordered qmatrix
all_snp_barplot_k5<- barplot(t(qmatrix2), col = c( "#FDDBC7","#4393C3","pink2", "#D1E5F0" , "#D6604D"),
                             border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients")+
  axis(1, at = seq_along(meta.data$site), labels = meta.data$site, las = 2, cex.axis = 0.7)

# Add site labels below the barplot
axis(1, at = seq_along(meta.data$site), labels = meta.data$site, las = 2, cex.axis = 0.7)


# Add sample names to qmatrix
rownames(qmatrix2) <- meta.data$sample

# Summarize qmatrix by site
qmatrix_by_site <- meta.data %>%
  select(sample, site2) %>%
  left_join(as.data.frame(qmatrix2) %>% rownames_to_column("sample"), by = "sample") 



# Aggregate proportions by site2
agg_data <- qmatrix_by_site %>%
  group_by(site2) %>%
  summarize(
    V1 = mean(V1),
    V2 = mean(V2),
    V3 = mean(V3),
    V4 = mean(V4)
  )

# Reshape data to long format
agg_data_long <- agg_data %>%
  pivot_longer(cols = V1:V4, names_to = "Population", values_to = "Proportion")




qmatrix_by_site_long <- qmatrix_by_site %>%
  pivot_longer(cols = starts_with("site2"), names_to = "population", values_to = "proportion")


# Plot pie charts for each site2
ggplot(agg_data_long, aes(x = "", y = Proportion, fill = Population)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  facet_wrap(~ site2, ncol = 3) + # Adjust the number of columns for layout
  scale_fill_manual(values = c("#D6604D", "#D1E5F0", "#4393C3", "#F4A582")) + # Customize colors
  labs(title = "Ancestry Proportions by Site2",
       fill = "Population") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "bottom"
  )



############################################
##################THINNED SNPS##############
############################################

#Genlight from previous file

gl_thinned  #138 samples, 4,171 variants


meta <- data.frame(site = c("BA", "BAD", "DD", "DS", "ED", "ES", "MBD", "MBS", "SQD", "SQS", "T"), 
                   site2 = c("BA", "BA","BS", "BS", "E", "E", "MB", "MB", "SQ", "SQ", "BS"),
                   depth = c("Shallow", "Deep", "Deep", "Shallow", "Deep", "Shallow", "Deep", "Shallow", "Deep", "Shallow", "Shallow")
)


meta


# Sample column names from vcf@gt
sample_names <- colnames(vcf_thinned@gt)

# Use gsub with a regex pattern to extract the first 1-3 letters, stopping at the first number
site <- gsub("^([A-Za-z]{1,3})\\d.*", "\\1", sample_names)[2:139]

meta2 <- data.frame(sample = sample_names[2:139], site = site)

meta.data <- left_join(meta, meta2, by = "site")

meta.data

gl_thinned

gl2faststructure(gl_thinned, outfile = "./gl_structure3.fstr",
                 outpath = './intermediate_files3/')


LEA::struct2geno("./intermediate_files3/gl_structure3.fstr", ploidy = 2, FORMAT = 2)
obj.snmf2 = snmf("./intermediate_files3/gl_structure3.fstr.geno", 1:10,
                entropy = TRUE, ploidy = 2,project = "new",  repetitions = 10, alpha = 100)
qmatrix2 = Q(obj.snmf2, K = 5, run = 10)

barplot(t(qmatrix2), col = c("#D1E5F0","#D6604D","#FDDBC7","#F4A582","#92C5DE","#4393C3"), border = NA, space = 0,
        xlab = "Individuals", ylab = "Admixture coefficients")

plot(obj.snmf2, pch = 19, cex = 1.5)


# Ensure meta.data has the correct mapping of samples and sites
meta.data <- left_join(meta2, meta, by = "site")

# Reorder meta.data by the desired site order
desired_order <- c("DD", "DS","T", "ED", "ES",  "SQD", "SQS", "BA", "BAD", "MBD", "MBS")
meta.data <- meta.data %>%
  mutate(site = factor(site, levels = desired_order)) %>%
  arrange(site)


# Reorder qmatrix rows to match the order of samples in meta.data
qmatrix3 <- qmatrix2[match(meta.data$sample, indNames(gl_thinned)), ]

# Create the barplot with the reordered qmatrix
barplot(t(qmatrix3), col = c("#D6604D","#D1E5F0", "#4393C3",  "#FDDBC7"),
        border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients")

# Optional: Add site labels below the barplot
axis(1, at = seq_along(meta.data$site), labels = meta.data$site, las = 2, cex.axis = 0.7)



#####Outlier SNPs only
gl_outliers <- gl_obj2[, outlier_ids]



gl_outliers  #138 samples, 1749 variants


meta <- data.frame(site = c("BA", "BAD", "DD", "DS", "ED", "ES", "MBD", "MBS", "SQD", "SQS", "T"), 
                   site2 = c("BA", "BA","BS", "BS", "E", "E", "MB", "MB", "SQ", "SQ", "BS"),
                   depth = c("Shallow", "Deep", "Deep", "Shallow", "Deep", "Shallow", "Deep", "Shallow", "Deep", "Shallow", "Shallow")
)


meta


# Sample column names from vcf@gt
sample_names <- gl_outliers@ind.names

# Use gsub with a regex pattern to extract the first 1-3 letters, stopping at the first number
site <- gsub("^([A-Za-z]{1,3})\\d.*", "\\1", sample_names)[2:139]

meta2 <- data.frame(sample = sample_names[2:139], site = site)

meta.data <- left_join(meta, meta2, by = "site")

meta.data

gl_outliers

gl2faststructure(gl_outliers, outfile = "./gl_structure.fstr",
                 outpath = './intermediate_files4/')


LEA::struct2geno("./intermediate_files4/gl_structure.fstr", ploidy = 2, FORMAT = 2)
obj.snmf3 = snmf("./intermediate_files4/gl_structure.fstr.geno", 1:10,
                 entropy = TRUE, ploidy = 2,project = "new",  repetitions = 10, alpha = 100)
qmatrix4 = Q(obj.snmf3, K = 5, run = 10)

barplot(t(qmatrix4), col = c("#D1E5F0","#D6604D","#FDDBC7","#F4A582","#92C5DE","#4393C3"), border = NA, space = 0,
        xlab = "Individuals", ylab = "Admixture coefficients")

plot(obj.snmf3, pch = 19, cex = 1.5)


# Ensure meta.data has the correct mapping of samples and sites
meta.data <- left_join(meta2, meta, by = "site")

# Reorder meta.data by the desired site order
desired_order <- c("DD", "DS","T", "ED", "ES",  "SQD", "SQS", "BA", "BAD", "MBD", "MBS")
meta.data <- meta.data %>%
  mutate(site = factor(site, levels = desired_order)) %>%
  arrange(site)


# Reorder qmatrix rows to match the order of samples in meta.data
qmatrix5 <- qmatrix4[match(meta.data$sample, indNames(gl_outliers)), ]

# Create the barplot with the reordered qmatrix
barplot(t(qmatrix5), col = c("#D6604D","#D1E5F0", "#4393C3",  "#FDDBC7", "pink"),
        border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients")

# Optional: Add site labels below the barplot
axis(1, at = seq_along(meta.data$site), labels = meta.data$site, las = 2, cex.axis = 0.7)





