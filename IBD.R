##Isolation by distance analysis
pairwise_fst # From the FST.R code file

## Isolation-by-distance using pairwise FST from hierfstat
## genetic distance = FST / (1 - FST)
## geographic distance = great-circle (Haversine) in km

# Packages
library(geosphere)  # distHaversine
library(vegan)      # mantel
library(ggplot2)

# ----------------------------
# 1) Site coordinates (as given)
# ----------------------------
coords <- data.frame(
  site = c("MBS","MBD","BA","BAD","SQS","SQD","ES","ED","T","DS","DD"),
  latitude  = c(24.628393,24.628393,27.129425,27.129425,30.513849,30.513849,31.855235,31.855235,48.827373,48.876831,48.877971),
  longitude = c(-112.158024,-112.158024,-114.288862,-114.288862,-116.067151,-116.067151,-116.655970,-116.655970,-125.197030,-125.092201,-125.087791)
)

# ----------------------------
# 2) Make sure pairwise_fst is a matrix with matching dimnames
# ----------------------------
fst_mat <- as.matrix(pairwise_fst$Fsts)

# Check that sites match between FST matrix and coords
sites <- intersect(rownames(fst_mat), coords$site)
if (length(sites) < 3) stop("FST matrix and coords have <3 overlapping sites. Check names.")

# Reorder both to the same site order
fst_mat <- fst_mat[sites, sites, drop = FALSE]
coords  <- coords[match(sites, coords$site), , drop = FALSE]

# ----------------------------
# 3) Genetic distance: FST/(1-FST)
# ----------------------------
# Optional: handle edge cases
# - if you have negative FST values, you can either keep them or clamp to 0 (common choice for IBD plots)
fst_use <- fst_mat
fst_use[is.na(fst_use)] <- NA
fst_use[fst_use >= 1] <- NA  # avoid division by zero / nonsense
# fst_use[fst_use < 0] <- 0  # uncomment if you want to clamp negatives to 0

gen_dist <- fst_use / (1 - fst_use)
diag(gen_dist) <- 0

# ----------------------------
# 4) Geographic distance matrix (km)
# ----------------------------
lonlat <- coords[, c("longitude", "latitude")]
geo_dist_m  <- geosphere::distm(lonlat, fun = geosphere::distHaversine)
geo_dist_km <- geo_dist_m / 1000
dimnames(geo_dist_km) <- list(coords$site, coords$site)

# ----------------------------
# 5) Mantel test (permutation-based)
# ----------------------------
# Mantel expects "dist" objects
mant <- vegan::mantel(as.dist(gen_dist), as.dist(geo_dist_km),
                      method = "pearson", permutations = 9999)
print(mant)

# ----------------------------
# 6) Simple IBD scatter + linear fit (non-independence warning!)
# ----------------------------
# Convert to pairwise table (lower triangle only)
to_pairs <- function(mat, value_name) {
  stopifnot(is.matrix(mat), nrow(mat) == ncol(mat))
  idx <- which(lower.tri(mat), arr.ind = TRUE)
  data.frame(
    site1 = rownames(mat)[idx[,1]],
    site2 = colnames(mat)[idx[,2]],
    value = mat[idx],
    stringsAsFactors = FALSE
  ) |>
    setNames(c("site1","site2", value_name))
}

df_g <- to_pairs(gen_dist, "fst_over_1minusfst")
df_d <- to_pairs(geo_dist_km, "dist_km")
ibd_df <- merge(df_g, df_d, by = c("site1","site2"))

# Remove NAs/infs if any
ibd_df <- subset(ibd_df, is.finite(fst_over_1minusfst) & is.finite(dist_km))

# Plot
ggplot(ibd_df, aes(x = dist_km, y = fst_over_1minusfst)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "darkgrey") +
  labs(x = "Geographic distance (km)",
       y = expression(F[ST] / (1 - F[ST])),
       title = "") +
  theme_classic()

# Optional: (naïve) lm summary (again: pairs are not independent)
summary(lm(fst_over_1minusfst ~ dist_km, data = ibd_df))

