# ==============================================================================
# File:         inla-spde.r
# Project:      Multi-sensor INLA-SPDE Bayesian Flood Mapping
# Author:       Francisco J. Lozada-Aguilar
# Affiliation:  Erasmus Mundus MSc Geospatial Technologies / UNGCS Collaboration
# Date:         2026-01-22
#
# Description:
#   This script implements a Bayesian spatial model for flood detection based on
#   a latent Gaussian field (LGF) constructed through the SPDE approach and
#   estimated using Integrated Nested Laplace Approximation (INLA)
#
#   The response variable is assumed to follow a Binomial/Bernoulli distribution,
#   while spatial dependence is captured by a continuously indexed Gaussian field
#   approximated as a Gaussian Markov Random Field (GMRF) over a triangulated mesh
#
#   Conceptually, the spatial field represents structured spatial variability
#   that cannot be explained solely by observed covariates (e.g., spectral indices,
#   topography), acting as a regularized latent process rather than a directly
#   observable physical quantity
#
# Statistical Model:
#   Let y_i ~ Binomial(n_i, p_i), where
#
#     logit(p_i) = η_i
#               = β_0 + Σ_k β_k x_{ik} + f(s_i)
#
#   where:
#     - β_k are fixed effects associated with observed covariates,
#     - f(s_i) is a latent spatial Gaussian field evaluated at location s_i,
#     - f(s) is defined as the solution to a stochastic partial differential
#       equation (SPDE), yielding a Matérn covariance structure
#
#   The continuous field f(s) is discretized via a finite element method (FEM)
#   on a triangular mesh, resulting in a sparse GMRF representation
#
# Computational Framework:
#   - SPDE formulation: Lindgren et al. (2011)
#   - Bayesian inference: INLA (Rue et al., 2009)
#   - Spatial domain discretization via mesh construction
#
# Inputs:
#   - response: Binary or binomial flood indicator (e.g., water / non-water)
#   - covariates: Remote sensing and geospatial predictors (NDVI, elevation, etc.)
#   - coordinates: Spatial locations in projected CRS
#   - mesh parameters: Controls spatial resolution and computational complexity
#
# Outputs:
#   - Posterior summaries of fixed effects (β)
#   - Posterior marginal distributions of SPDE hyperparameters (range, variance)
#   - Spatial predictions of flood probability
#   - Latent field estimates for spatial uncertainty analysis
#
# Methodological Notes:
#   - The mesh design determines the effective spatial scale of inference and
#     should reflect both sensor resolution and process smoothness
#   - Priors on SPDE hyperparameters encode assumptions about spatial correlation
#     length and marginal variance
#   - The latent field captures residual spatial structure and should not be
#     interpreted as direct physical causation
#
# Limitations:
#   - Assumes conditional independence of observations given the latent field
#   - Stationarity and isotropy are imposed through the Matérn covariance
#   - Model performance is sensitive to prior choice and mesh configuration
#
# Reproducibility:
#   - R version: <R version>
#   - Required packages: INLA, sp, sf, raster, terra
#   - Random seeds are fixed where stochastic components are involved
#
# References:
#   - Lindgren, F., Rue, H., & Lindström, J. (2011).
#     An explicit link between Gaussian fields and Gaussian Markov random fields
#     Journal of the Royal Statistical Society: Series B.
#   - Rue, H., Martino, S., & Chopin, N. (2009)
#     Approximate Bayesian inference for latent Gaussian models using INLA
#
# License:
#   MIT License
#   Copyright (c) 2026 Francisco J. Lozada-Aguilar
# ============================================================================== 


# ==============================================================================
# 0. REQUIRED LIBRERIES
# ==============================================================================
library(sp)
library(terra)
library(sf)
library(raster)
library(INLA)
library(ggplot2)
library(viridis)
library(pROC)

# ==============================================================================
# 1. AUXILIARY FUNCTIONS
# ==============================================================================

extract_dates <- function(path)
{
  files <- list.files(path, pattern = "\\.shp$|\\.SHP$", full.names = FALSE)
  pattern <- "(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dic)[0-9]{4}"
  dates <- regmatches(files, gregexpr(pattern, files, ignore.case = TRUE))
  return(tolower(unique(unlist(dates))))
}

align_raster <- function(raster_input, raster_master)
{
  r_in <- rast(raster_input)
  if (!compareGeom(r_in, raster_master, stopOnError = FALSE)) {
    if (crs(r_in) != crs(raster_master)) {
      r_in <- project(r_in, raster_master)
    }
    r_in <- resample(r_in, raster_master, method = "bilinear")
  }
  return(r_in)
}

normalize_raster_terra <- function(r_input) 
{
  mm <- minmax(r_input)
  min_val <- mm[1, ]
  max_val <- mm[2, ]
  if (max_val == min_val) return(r_input * 0)
  r_norm <- (r_input - min_val) / (max_val - min_val)
  return(r_norm)
}

# ==============================================================================
# 2. EXECUTION PARAMETERS AND ROUTES
# ==============================================================================

AOI  <- 6
path      <- paste0("data/trainingdata/AOI", AOI, "/")
pathTest  <- paste0("data/testdata/AOI", AOI, "/")
pathResults <-paste0("results/AOI", AOI, "/")

dates       <- extract_dates(path)
datesTest   <- extract_dates(pathTest)
master_ref  <- rast(paste0(path, "s1vv-", dates[1], ".tif"))

# ==============================================================================
# 3. STATIC RASTERS
# ==============================================================================
dem   <- align_raster(paste0(path, "dem.tif"),        master_ref)[[1]]
slope <- align_raster(paste0(path, "slope.tif"),      master_ref)
hand  <- align_raster(paste0(path, "hand.tif"),       master_ref)
curv  <- align_raster(paste0(path, "curvature.tif"),  master_ref)
twi   <- align_raster(paste0(path, "twi.tif"),        master_ref)

names(dem)    <- "dem"
names(slope)  <- "slope"
names(hand)   <- "hand"
names(curv)   <- "curv"
names(twi)    <- "twi"

# ==============================================================================
# 4. DYNAMIC RASTERS + POINT PATTERNS
# ==============================================================================
all_points_sf      <- list()
dynamic_stack_list <- list()

cat("Processing temporal fusion...\n")

for (date in dates) 
{
  pp_temp <- st_read(paste0(path, "pp-", date, ".shp"), quiet = TRUE)
  pp_temp <- st_transform(pp_temp, crs(master_ref))
  all_points_sf[[date]] <- pp_temp
  
  s1_vv_t <- align_raster(paste0(path, "s1vv-",  date, ".tif"), master_ref)
  s1_vh_t <- align_raster(paste0(path, "s1vh-",  date, ".tif"), master_ref)
  ndvi_t  <- align_raster(paste0(path, "ndvi-",  date, ".tif"), master_ref)
  mndwi_t <- align_raster(paste0(path, "mndwi-", date, ".tif"), master_ref)
  awei_t  <- align_raster(paste0(path, "awei-",  date, ".tif"), master_ref) 
  
  r_stack <- c(s1_vv_t, s1_vh_t, ndvi_t, mndwi_t, awei_t)
  names(r_stack) <- c("vv", "vh", "ndvi", "mndwi", "awei")
  dynamic_stack_list[[date]] <- r_stack
}

merged_sf <- do.call(rbind, all_points_sf)

# ==============================================================================
# 5. TEMPORAL AGGREGATION AND STANDARDIZATION
# ==============================================================================
big_stack <- rast(dynamic_stack_list)
num_vars  <- 5
num_dates <- nlyr(big_stack) / num_vars
indices   <- rep(1:num_vars, times = num_dates)

dynamic_mean <- tapp(big_stack, index = indices, fun = mean, na.rm = TRUE)
names(dynamic_mean) <- c("vv", "vh", "ndvi", "mndwi", "awei")

raw_covariates <- c(dem, slope, hand, curv, twi, dynamic_mean)

cat("Scaling rasters (Z-score)...\n")
final_covariates_terra <- scale(raw_covariates)
names(final_covariates_terra) <- names(raw_covariates)

# ==============================================================================
# 6. MESH AND SPDE SETTING
# ==============================================================================
locs_obs <- st_coordinates(merged_sf)

# Calculate diagonal distance to gauge scale
diag_dist <- sqrt(diff(ext(final_covariates_terra)[1:2])^2 + diff(ext(final_covariates_terra)[3:4])^2)

# Mesh settings 
# Finer resolution inside
max_edge_inner <- diag_dist / 150
# Coarser outside
max_edge_outer <- diag_dist / 20

mesh <- inla.mesh.2d(
  loc = locs_obs,
  max.edge = c(max_edge_inner, max_edge_outer),
  cutoff   = max_edge_inner / 2, # Avoid points too close
  offset   = c(max_edge_inner, max_edge_outer)
)

cat("Mesh Nodes:", mesh$n, "\n")

# SPDE Definition
prior_range_val <- diag_dist * 0.2 

spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  # P(range < prior_val) = 0.01
  prior.range = c(prior_range_val, 0.01),
  # P(sigma > 1) = 0.01
  prior.sigma = c(1, 0.01)                
)

# ==============================================================================
# 7. DATASET PREPARATION (BINARY RESPONSE + STACK INLA)
# ==============================================================================

cat("Preparing data for INLA...\n")

merged_v <- terra::vect(merged_sf)

if (!terra::same.crs(merged_v, final_covariates_terra))
{
  merged_v <- terra::project(merged_v, terra::crs(final_covariates_terra))
}

r_y <- terra::rasterize(
  merged_v,
  final_covariates_terra[[1]],
  field = 1,
  fun = "max",
  background = 0
)
names(r_y) <- "flood_obs"

# Full Stack  (response + covariates) ---
full_stack <- c(r_y, final_covariates_terra)

# DataFrame (xy + variables), cleaning NA's ---
df_pixels_raw <- terra::as.data.frame(full_stack, xy = TRUE, na.rm = TRUE)

n_data <- nrow(df_pixels_raw)
cat("Number of valid pixels:", n_data, "\n")
cat("Verification flood_obs (must be 0/1):\n")
print(table(df_pixels_raw$flood_obs))

req_vars <- c("dem","slope","hand","curv","twi","vv", "vh","ndvi","mndwi","awei")
missing_vars <- setdiff(req_vars, names(df_pixels_raw))
if (length(missing_vars) > 0) 
{
  stop("Missing covariates in df_pixels_raw: ", paste(missing_vars, collapse = ", "))
}

df_pixels_raw[req_vars] <- scale(df_pixels_raw[req_vars])

# ==============================================================================
# 8. STACK BUILDING
# ==============================================================================

cat("CBuilding INLA matrices...\n")

loc_pixels <- as.matrix(df_pixels_raw[, c("x", "y")])
A_matrix_raw <- inla.spde.make.A(mesh, loc = loc_pixels)

idx_spatial <- inla.spde.make.index("s", n.spde = spde$n.spde)


stk_pix <- inla.stack(
  data = list(
    y = df_pixels_raw$flood_obs,
    Ntrials = rep(1, n_data)   # Bernoulli
  ),
  A = list(1, A_matrix_raw),
  effects = list(
    data.frame(
      Intercept = 1,
      dem  = df_pixels_raw$dem,
      slope = df_pixels_raw$slope,
      hand = df_pixels_raw$hand,
      curv  = df_pixels_raw$curv,
      twi  = df_pixels_raw$twi,
      vv    = df_pixels_raw$vv,
      vh    = df_pixels_raw$vh,
      ndvi  = df_pixels_raw$ndvi,
      mndwi = df_pixels_raw$mndwi,
      awei = df_pixels_raw$awei
    ),
    idx_spatial
  ),
  tag = "pixels"
)

df_pixels <- inla.stack.data(stk_pix)
A_pixels  <- inla.stack.A(stk_pix)

cat("Data ready\n")

# ==============================================================================
# 9. IMPLEMENTATION OF THE INLA MODEL
# ==============================================================================

# Formula
formula <- y ~ 0 + Intercept + dem + slope + hand  + curv + twi + vv + vh + ndvi + mndwi + awei + f(s, model = spde)

cat("Running INLA (Binomial)...\n")

result <- inla(
  formula,
  family = "binomial",
  data = df_pixels,
  Ntrials = df_pixels$Ntrials,
  control.predictor = list(A = A_pixels, compute = TRUE),
  
  control.inla = list(
    strategy = "adaptive",        # Smart optimization step sizing
    int.strategy = "eb",          # Empirical Bayes
    control.vb = list(enable = FALSE), # Disables the Variational Bayes correction that caused the crash
    tolerance = 1e-5              # Tolerance slightly to aid convergence
  ),
  # --------------------------------
  
  control.compute = list(waic = TRUE, cpo = TRUE, config = TRUE),
  verbose = TRUE,
  safe = TRUE
)

summary(result)

# ==============================================================================
# 10. RESULTS AND DATA VISUALIZATION
# ==============================================================================

cat("\nSummary of Fixed Effects:\n")
print(result$summary.fixed)

cat("Generating a susceptibility map...\n")

idx_pixels <- inla.stack.index(stk_pix, tag = "pixels")$data
prob_flood <- result$summary.fitted.values$mean[idx_pixels]
df_plot <- df_pixels_raw[, c("x", "y")]
df_plot$probability <- prob_flood
r_final <- rast(df_plot, type="xyz", crs=crs(final_covariates_terra))

# Export raster
writeRaster(r_final, paste0(pathResults, "Flood_Susceptibility_Binomial.tif"), overwrite=TRUE)

# Visualization with ggplot
ggplot(df_plot, aes(x=x, y=y, fill=probability)) +
  geom_raster() +
  scale_fill_viridis_c(option="turbo", name="Flood\nProb.", limits=c(0,1)) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Flood Susceptibility Map",
       subtitle = "Spatial Logistic Regression (INLA-SPDE)",
       x = "Longitude", y = "Latitude")


# SPATIAL LATENT FIELD

gproj <- inla.mesh.projector(mesh, dims = c(300, 300))
valores_campo_latente <- result$summary.random$s$mean
spatial_mean_proj <- inla.mesh.project(gproj, valores_campo_latente)
r_spatial <- rast(list(x = gproj$x, 
                       y = gproj$y, 
                       z = spatial_mean_proj)) 

crs(r_spatial) <- crs(final_covariates_terra)
plot(r_spatial, main="Spatial Latent Field", col=viridis::viridis(100))
writeRaster(r_final, paste0(pathResults, "latent_field.tif"), overwrite=TRUE)

# COVARIATES EFFECT 
fixed_effects <- result$summary.fixed
fixed_effects$variable <- rownames(fixed_effects)
fixed_effects <- fixed_effects[fixed_effects$variable != "Intercept", ]

ggplot(fixed_effects, aes(x = mean, y = reorder(variable, mean))) +
  geom_point(size = 3, color = "blue") +
  geom_errorbarh(aes(xmin = `0.025quant`, xmax = `0.975quant`), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Covariates effect",
       x = "Size effect (Betha)", y = "Covariable") +
  theme_bw()

# ROC AND AUC VALIDATION
roc_obj <- roc(df_pixels_raw$flood_obs, df_plot$probability)

print(paste("AUC:", auc(roc_obj)))
plot(roc_obj, main="ROC Curve - Flood susceptibility")

# ==============================================================================
# 11. TESTING STAGE
# ==============================================================================

cat("\n====================\n")
cat("TESTING...\n")
cat("====================\n")

get_train_scale_stats <- function(raw_covs, vars_keep) 
{
  stats <- data.frame(
    var  = names(raw_covs),
    mean = as.numeric(terra::global(raw_covs, "mean", na.rm = TRUE)[,1]),
    sd   = as.numeric(terra::global(raw_covs, "sd",   na.rm = TRUE)[,1])
  )
  stats <- stats[stats$var %in% vars_keep, ]
  stats$sd[is.na(stats$sd) | stats$sd == 0] <- 1
  rownames(stats) <- stats$var
  stats
}

vars_model <- c("dem","slope","hand","curv","twi","vv", "vh","ndvi","mndwi","awei")
train_stats <- get_train_scale_stats(raw_covariates, vars_model)
print(train_stats)

# Date prediction function (with chunking)

predict_flood_date <- function(date_str, path_data,
                               model_inla, mesh, train_stats,
                               master_ref, dem, slope, hand, curv, twi,
                               chunk_size = 200000) 
{
  
  message(paste0("\n--- TEST date: ", date_str, " ---"))
  
  betas <- setNames(model_inla$summary.fixed$mean, rownames(model_inla$summary.fixed))
  
  if (!("Intercept" %in% names(betas)) && ("(Intercept)" %in% names(betas))) {
    betas["Intercept"] <- betas["(Intercept)"]
  }
  if (!("Intercept" %in% names(betas)))
  {
    stop("Intercept not found in summary.fixed. Available names: ",
         paste(names(betas), collapse=", "))
  }
  
  s_mean <- model_inla$summary.random$s$mean
  if (is.null(s_mean)) stop("Model_inla$summary.random$s$mean not found")
  
  f_vv    <- paste0(path_data, "s1vv-",  date_str, ".tif")
  f_vh    <- paste0(path_data, "s1vh-",  date_str, ".tif")
  f_ndvi  <- paste0(path_data, "ndvi-",  date_str, ".tif")
  f_mndwi <- paste0(path_data, "mndwi-", date_str, ".tif")
  f_awei  <- paste0(path_data, "awei-",  date_str, ".tif")
  
  for (ff in c(f_vv, f_vh, f_ndvi, f_mndwi, f_awei)) 
  {
    if (!file.exists(ff)) stop("File not found: ", ff)
  }
  
  r_vv    <- align_raster(f_vv,   master_ref)
  r_vh    <- align_raster(f_vh,   master_ref)
  r_ndvi  <- align_raster(f_ndvi, master_ref)
  r_mndwi <- align_raster(f_mndwi,master_ref)
  r_awei  <- align_raster(f_awei, master_ref)
  
  names(r_vv) <- "vv"; names(r_vh) <- "vh"; names(r_ndvi) <- "ndvi"; names(r_mndwi) <- "mndwi"; names(r_awei) <- "awei"
  
  vars_model <- c("dem","slope","hand","curv","twi","vv","vh","ndvi","mndwi","awei")
  
  stack_test <- c(dem, slope, hand, curv, twi, r_vv, r_vh, r_ndvi, r_mndwi, r_awei)
  names(stack_test) <- vars_model
  
  # Pixel dataframe
  df_test <- terra::as.data.frame(stack_test, xy = TRUE, na.rm = TRUE)
  if (nrow(df_test) == 0) stop("df_test quedó vacío (todo NA).")
  n <- nrow(df_test)
  cat("Píxeles válidos en test:", n, "\n")
  
  for (v in vars_model) 
  {
    mu <- train_stats[v, "mean"]
    sd <- train_stats[v, "sd"]
    if (is.na(sd) || sd == 0) sd <- 1
    df_test[[v]] <- (df_test[[v]] - mu) / sd
  }
  
  # etha = Intercept + Xb + A*s
  eta <- rep(betas["Intercept"], n)
  
  covs_in_model <- intersect(vars_model, names(betas))
  covs_in_model <- setdiff(covs_in_model, c("Intercept","(Intercept)"))
  
  if (length(covs_in_model) > 0) {
    X <- as.matrix(df_test[, covs_in_model, drop = FALSE])
    b <- betas[covs_in_model]
    eta <- eta + as.vector(X %*% b)
  } else {
    warning("Covariates not found (other than Intercept) in betas. Beta names: ",
            paste(names(betas), collapse=", "))
  }
  
  coords <- as.matrix(df_test[, c("x","y")])
  idx <- seq_len(n)
  chunks <- split(idx, ceiling(idx / chunk_size))
  cat("Chunks:", length(chunks), " (chunk_size=", chunk_size, ")\n", sep="")
  
  spatial_term <- numeric(n)
  for (k in seq_along(chunks)) {
    ii <- chunks[[k]]
    A_pred <- inla.spde.make.A(mesh, loc = coords[ii, , drop = FALSE])
    spatial_term[ii] <- as.vector(A_pred %*% s_mean)
    if (k %% 5 == 0) cat("  ...chunk", k, "de", length(chunks), "\n")
  }
  
  eta <- eta + spatial_term
  p <- plogis(eta)
  df_out <- df_test[, c("x","y")]
  df_out$probability <- p
  
  r_prob <- terra::rast(df_out, type = "xyz", crs = terra::crs(master_ref))
  r_prob <- terra::resample(r_prob, master_ref, method = "bilinear")
  r_prob <- terra::mask(r_prob, dem)
  
  names(r_prob) <- paste0("prob_", date_str)
  r_prob
}

test_results <- list()

for (d in datesTest) {
  
  p_map <- predict_flood_date(
    date_str   = d,
    path_data  = pathTest,
    model_inla = result,
    mesh       = mesh,
    train_stats = train_stats,
    master_ref = master_ref,
    dem = dem, slope = slope, hand = hand, curv = curv, twi= twi,
    chunk_size = 200000
  )
  
  out_file <- paste0(pathResults, "Pred_Test_", d, ".tif")
  terra::writeRaster(p_map, out_file, overwrite = TRUE)
  
  test_results[[d]] <- p_map
  
  plot(p_map, main = paste("Flood probability -", d),
       range = c(0,1), col = terra::map.pal("viridis"))
}