# EiseniaPopGen

Code used for population genomic and seascape genomic analyses in:

**Starko, S.** et al. (2026) *Environmental Drivers of Genetic Structure and Local Adaptation in a Marine Foundation Species* (Ecology & Evolution).

This repo contains R scripts used to generate the population structure, diversity, connectivity, and genotype–environment association results for *Eisenia arborea* ddRAD SNP data.

> Note: The scripts were written as an analysis “notebook in .R files” and assume some objects are created in a prior script / interactive session. 
## Contents

- `vcf_filtered2.vcf`  
  Filtered SNP dataset used for the majority of downstream analyses (9315 SNPs / 138 individuals in the final dataset in the paper).
- R scripts (high-level):
  - SNP thinning / prep: `vcf_thinning.R`, `SFS_prep3.R`
  - Population structure: `snmf.R`, `PCA_no_outliers.R`
  - Differentiation & connectivity: `FST.R`, `IBD.R`, `migration.R`
  - Genetic diversity / inbreeding: `Inbreeding_stats.R`
  - Ne / LDNe: `Ne-estimator.R`
  - Genotype–environment association: `RDA.R`, `lfmm_depth.R`
  - Error rates / replicates: `ErrorRates.R`
  - Plotting external demographic outputs: `PlotStairwayPlotResults.R`
  - Utility: `CreateGenpopfile.R`, `Pie_charts_forMainFig.R`

## Requirements

Developed and run in **R 4.3.2**.

Main R packages used across scripts include (not exhaustive):
- `vcfR`, `adegenet`, `dartR.base` / `dartR.popgen`
- `SNPRelate`, `LEA` / `snmf`
- `hierfstat`, `StAMPP`, `diveRsity`, `divMigrateR`
- `vegan`, `lfmm`, `qvalue`
- `geosphere`, `ggplot2`, `dplyr`, `reshape2`

### External tools (upstream / demographic workflows)
Some steps referenced in the paper (e.g., STACKS, vcftools, BayeScan, StairwayPlot2, EPOS, easySFS) are not run directly from a single “pipeline” script here. Where outputs from these tools are needed for plotting, the scripts assume those result files exist in the expected folders.
