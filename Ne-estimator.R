library(dartR.popgen)

gl_obj2

gl_thinned <- vcfR2genlight(vcf_thinned)
gl_thinned@pop <- as.factor(extracted_sites[2:139])
gl_thinned2 <- gl_thinned[pop(gl_thinned) == "BA"]

run2 <- gl.LDNe(gl_obj2, outfile = "NeEstimate_all.txt", critical = c(0.03, 0.05), mating = "random", pairing = "separate")

extracted_sites




# install.packages("ggplot2")  # if not already installed
library(ggplot2)

# 1. Create your data frame (unchanged order in the vector):
df <- data.frame(
  Population = c("BA", "DS", "ED", "ES", "MBS", "SQD", "SQS"),
  Ne         = c(60.7, 14.2, 87.6, 16.3, 13.3, 3.3, 13.8),
  lower_ci   = c(57.5, 13.7, 81.2, 16.0, 12.8, 3.3, 13.3),
  upper_ci   = c(64.3, 14.7, 95.1, 16.6, 13.9, 3.6, 14.4)
)

# 2. Plot using ggplot with a custom x-axis order
ggplot(df, aes(x = Population, y = Ne)) +
  geom_col(fill = "skyblue") +
  geom_errorbar(
    aes(ymin = lower_ci, ymax = upper_ci),
    width = 0.2,
    color = "black"
  ) +
  # Specify the order of categories on the x-axis:
  scale_x_discrete(limits = c("MBS", "BA", "SQD", "SQS", "ES", "ED", "DS")) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "",
    x = "Population",
    y = "Effective Population Size (Ne)"
  )


# Make sure Population is a factor with the specified levels
df$Population <- factor(df$Population, levels = c("MBS", "BA", "SQD", "SQS", "ES", "ED", "DS"))

plot_Ne <- ggplot(df, aes(x = Population, y = Ne)) +
  geom_col(fill = "grey") +
  geom_errorbar(
    aes(ymin = lower_ci, ymax = upper_ci),
    width = 0.2,
    color = "black"
  ) +
  scale_x_discrete(limits = c("MBS", "MBD", "BA", "BAD", "SQS", "SQD", "ES", "ED", "T", "DS", "DD")) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(
    title = "Effective population size (Ne)",
    x = "Population",
    y = "Effective Population Size (Ne)"
  )





###StairwayPlot implemented in dartRverse

gl_thinned
gl_thinned2 <- gl_thinned[pop(gl_thinned) == "ED"]

#number of individuals
nInd(gl_thinned2)

#number of loci
nLoc(gl_thinned2)

gl_imp <- gl.impute(gl_thinned2)


#create 
sfscon <- gl.sfs(gl_imp, folded = TRUE)
gl.download.binary(software="stairway2",os="windows")

#L <- 4171
#mu = 1e-8
#takes about 5 minutes so not run here
#DS
#Ne_sw <- gl.run.stairway2(gl_imp, stairway2.path = file.path(tempdir(),"stairway2"), mu = mu, gentime = 2, run=TRUE, nreps = 200, parallel=2, L=L, minbinsize =1, plot.theme = theme_classic())
#Ne_sw_BA <- gl.run.stairway2(gl_imp, stairway2.path = file.path(tempdir(),"stairway2"), mu = mu, gentime = 2, run=TRUE, nreps = 200, parallel=2, L=L, minbinsize =1, plot.theme = theme_classic())
#Ne_sw_MBS <- gl.run.stairway2(gl_imp, stairway2.path = file.path(tempdir(),"stairway2"), mu = mu, gentime = 2, run=TRUE, nreps = 200, parallel=2, L=L, minbinsize =1, plot.theme = theme_classic())
#Ne_sw_SQS <- gl.run.stairway2(gl_imp, stairway2.path = file.path(tempdir(),"stairway2"), mu = mu, gentime = 2, run=TRUE, nreps = 200, parallel=2, L=L, minbinsize =1, plot.theme = theme_classic())
#Ne_sw_SQD <- gl.run.stairway2(gl_imp, stairway2.path = file.path(tempdir(),"stairway2"), mu = mu, gentime = 2, run=TRUE, nreps = 200, parallel=2, L=L, minbinsize =1, plot.theme = theme_classic())
#Ne_sw_ES <- gl.run.stairway2(gl_imp, stairway2.path = file.path(tempdir(),"stairway2"), mu = mu, gentime = 2, run=TRUE, nreps = 200, parallel=2, L=L, minbinsize =1, plot.theme = theme_classic())
#Ne_sw_ED <- gl.run.stairway2(gl_imp, stairway2.path = file.path(tempdir(),"stairway2"), mu = mu, gentime = 2, run=TRUE, nreps = 200, parallel=2, L=L, minbinsize =1, plot.theme = theme_classic())



mu = 1e-8
L =550000

vcf_thinned


# Get the SNP IDs that are flagged as outliers
outlier_ids <- snp_summary4$snpID[snp_summary4$outlier]

# Identify the row indices in the vcfR object that match these IDs.
# Adjust the column name ("ID") if your vcfR object stores locus names in a different column.
#rows_to_keep <- which(!vcf_thinned@fix[,"ID"] %in% outlier_ids)


####EPOS
gl_no_outliers <- gl_thinned[, !(locNames(gl_thinned) %in% outlier_ids)]
gl_thinned2 <- gl_no_outliers[pop(gl_no_outliers) == "BA"]

gl.download.binary("epos",os="windows")

epos <- gl.run.epos(gl_thinned2, epos.path = file.path(tempdir(),"epos"), L=4171, u = 1.25e-8, boot = 100, minbinsize = 1)
epos$history
BA_epos <- ggplot(epos$history, aes(x=generation, y=median))+geom_line()+
geom_ribbon(aes(ymin=low, ymax=high), alpha=0.2)+theme_classic()+ylim(c(0,6e6))+xlim(c(0, 6e6))+ylab("Effective population size (Ne)")


gl_thinned2 <- gl_thinned[pop(gl_thinned) == "SQS"]
gl.download.binary("epos",os="windows")
epos <- gl.run.epos(gl_thinned2, epos.path = file.path(tempdir(),"epos"), L=4171, u = 1.25e-8, boot = 100, minbinsize = 1)
epos$history
SQS_epos <- ggplot(epos$history, aes(x=generation, y=median))+geom_line()+
  geom_ribbon(aes(ymin=low, ymax=high), alpha=0.2)+theme_classic()+ylim(c(0,6e6))+xlim(c(0, 6e6))+ylab("Effective population size (Ne)")


gl_thinned2 <- gl_thinned[pop(gl_thinned) == "SQD"]
gl.download.binary("epos",os="windows")
epos <- gl.run.epos(gl_thinned2, epos.path = file.path(tempdir(),"epos"), L=4171, u = 1.25e-8, boot = 100, minbinsize = 1)
epos$history
SQD_epos <- ggplot(epos$history, aes(x=generation, y=median))+geom_line()+
  geom_ribbon(aes(ymin=low, ymax=high), alpha=0.2)+theme_classic()+ylim(c(0,6e6))+xlim(c(0, 6e6))+ylab("Effective population size (Ne)")

gl_thinned2 <- gl_thinned[pop(gl_thinned) == "ES"]
gl.download.binary("epos",os="windows")
epos <- gl.run.epos(gl_thinned2, epos.path = file.path(tempdir(),"epos"), L=4171, u = 1.25e-8, boot = 100, minbinsize = 1)
epos$history
ES_epos <- ggplot(epos$history, aes(x=generation, y=median))+geom_line()+
  geom_ribbon(aes(ymin=low, ymax=high), alpha=0.2)+theme_classic()+ylim(c(0,6e6))+xlim(c(0, 6e6))+ylab("Effective population size (Ne)")

gl_thinned2 <- gl_thinned[pop(gl_thinned) == "ED"]
gl.download.binary("epos",os="windows")
epos <- gl.run.epos(gl_thinned2, epos.path = file.path(tempdir(),"epos"), L=4171, u = 1.25e-8, boot = 100, minbinsize = 1)
epos$history
ED_epos <- ggplot(epos$history, aes(x=generation, y=median))+geom_line()+
  geom_ribbon(aes(ymin=low, ymax=high), alpha=0.2)+theme_classic()+ylim(c(0,6e6))+xlim(c(0, 6e6))+ylab("Effective population size (Ne)")

gl_thinned2 <- gl_thinned[pop(gl_thinned) == "DS"]
gl.download.binary("epos",os="windows")
epos <- gl.run.epos(gl_thinned2, epos.path = file.path(tempdir(),"epos"), L=4171, u = 1.25e-8, boot = 100, minbinsize = 1)
epos$history
DS_epos <- ggplot(epos$history, aes(x=generation, y=median))+geom_line()+
  geom_ribbon(aes(ymin=low, ymax=high), alpha=0.2)+theme_classic()+ylim(c(0,6e6))+xlim(c(0, 6e6))+ylab("Effective population size (Ne)")

gl_thinned2 <- gl_thinned[pop(gl_thinned) == "MBS"]
gl.download.binary("epos",os="windows")
epos <- gl.run.epos(gl_thinned2, epos.path = file.path(tempdir(),"epos"), L=4171, u = 1.25e-8, boot = 100, minbinsize = 1)
epos$history
MBS_epos <- ggplot(epos$history, aes(x=generation, y=median))+geom_line()+
  geom_ribbon(aes(ymin=low, ymax=high), alpha=0.2)+theme_classic()+ylim(c(0,6e6))+xlim(c(0, 6e6))+ylab("Effective population size (Ne)")

plot_grid(MBS_epos, BA_epos, SQD_epos, SQS_epos, ED_epos, ES_epos, DS_epos)

