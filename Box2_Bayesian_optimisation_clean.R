# ==============================================================================
# Box 2: Bayesian optimisation of the hybrid SDM-metapopulation model
# ==============================================================================
# Purpose
#   Reproduce the Box 2 simulation, Bayesian optimisation, validation, and
#   publication figures used to examine parameter recovery and equifinality in
#   the hybrid SDM-metapopulation model.
#
# Workflow
#   1. Generate a reference occupancy surface from known parameters.
#   2. Run 20 independent Bayesian optimisation runs against the complete
#      reference occupancy surface.
#   3. Select the parameter combination with the highest Warren's I.
#   4. Re-run the hybrid model 50 times using that parameterisation.
#   5. Compare the mean predicted surface with the reference surface.
#   6. Export the occupancy comparison (Figure 1) and parameter estimates
#      across optimisation runs (Figure 2).
#
# Notes
#   - Candidate parameterisations are evaluated against the complete reference
#     occupancy surface rather than a sampled subset.
#   - The script assumes that `landcover` and `veghgt` are available after
#     loading PopScape (or have been created before running the analysis).
#   - The C++ source file is expected at src/PopScapecpp.cpp by default. Set the
#     environment variable POPSCAPE_CPP to use a different location.
#   - Output paths are relative to the project root so that the script is
#     portable across machines.
# ==============================================================================


# ==============================================================================
# 1. PACKAGES AND PROJECT SETTINGS
# ==============================================================================

library(PopScape)

# Reproducible random-number sequence for simulation and optimisation.
set.seed(123)

# Relative output directories.
results_dir <- file.path("outputs", "results")
figures_dir <- file.path("outputs", "figures")
cache_dir <- file.path("cache", "natural_earth")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# Compile the C++ functions required by the metapopulation model.
cpp_file <- Sys.getenv(
  "POPSCAPE_CPP",
  unset = file.path("src", "PopScapecpp.cpp")
)

if (!file.exists(cpp_file)) {
  stop(
    "C++ source file not found: ", cpp_file, "\n",
    "Place PopScapecpp.cpp in src/ or set the POPSCAPE_CPP environment variable."
  )
}

Rcpp::sourceCpp(cpp_file)

# Check that the habitat data required below are available.
required_objects <- c("landcover", "veghgt")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]

if (length(missing_objects) > 0) {
  stop(
    "Required input object(s) not found: ",
    paste(missing_objects, collapse = ", "),
    "."
  )
}


# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

# Generate the reference occupancy surface from a known parameterisation.
generate_reference_surface <- function(
    habclass1,
    habclass2,
    habcont,
    msk,
    p1,
    p2,
    p3,
    p4,
    mu,
    x,
    alpha,
    ygamma,
    aos,
    timesteps = 500,
    minden = 0.01,
    maxden = 20,
    asprob = TRUE,
    n_runs = 50,
    verbose = TRUE
) {
  # Habitat suitability from the SDM component.
  linear_predictor <-
    p1 * habclass1 +
    p2 * habclass2 +
    p3 * log(habcont) +
    p4

  hsuit <- terra::mask(
    1 / (1 + exp(-linear_predictor)),
    msk
  )

  if (all(!is.finite(terra::values(hsuit)))) {
    stop("Reference habitat-suitability surface contains no finite values.")
  }

  occupancy_runs <- vector("list", n_runs)

  for (i in seq_len(n_runs)) {
    if (verbose) {
      message("Reference simulation ", i, " of ", n_runs)
    }

    occupancy_runs[[i]] <- MetaPopSim(
      hsuit,
      mu = mu,
      x = x,
      alpha = alpha,
      ygamma = ygamma,
      aos = aos,
      timesteps = timesteps,
      minden = minden,
      maxden = maxden,
      asprob = asprob
    )
  }

  occupancy_stack <- terra::rast(occupancy_runs)
  mean(occupancy_stack, na.rm = TRUE)
}


# Run one realisation of the hybrid SDM-metapopulation model.
OptimMetaPopHybrid <- function(
    habclass1,
    habclass2,
    habcont,
    msk,
    p1,
    p2,
    p3,
    p4,
    mu,
    x,
    alpha,
    ygamma,
    aos,
    timesteps = 500,
    minden = 0.01,
    maxden = 20
) {
  # Habitat suitability from the SDM component.
  linear_predictor <-
    p1 * habclass1 +
    p2 * habclass2 +
    p3 * log(habcont) +
    p4

  hsuit <- terra::mask(
    1 / (1 + exp(-linear_predictor)),
    msk
  )

  if (all(!is.finite(terra::values(hsuit)))) {
    return(NULL)
  }

  # Remove unsuitable habitat and delineate contiguous patches.
  hsuit[hsuit == 0] <- NA
  habitat <- hsuit * 0 + 1

  patch <- terra::patches(habitat, directions = 8)
  patch_matrix <- terra::as.matrix(patch, wide = TRUE)
  patch_matrix <- renumberPatchIds(patch_matrix)
  patch_matrix[patch_matrix == 0] <- NA

  # Convert habitat suitability to relative population density.
  rden <- log(hsuit / (1 - hsuit))^aos
  rden[is.na(rden)] <- minden
  rden[rden < minden] <- minden
  rden[rden > maxden] <- maxden
  rden_matrix <- terra::as.matrix(rden, wide = TRUE)

  # Patch-to-patch distances in kilometres.
  dij <- patchdistcpp(patch_matrix, rden_matrix) * terra::res(patch)[1] / 1000

  # Run the stochastic metapopulation model.
  Oij <- runmetapopmodel(
    patch_matrix,
    rden_matrix,
    terra::res(patch)[1],
    dij,
    mu,
    x,
    alpha,
    ygamma,
    timesteps = timesteps
  )

  # Convert patch-level model output back to a raster occupancy surface.
  patch_area <-
    (calculatePatchSizes(patch_matrix) * terra::res(patch)[1]^2) /
    (100 * 100)

  area_raster <- .rast(assignPatchValues(patch_matrix, patch_area), patch)

  extinction_probability <- mu / (area_raster * rden)^x
  extinction_probability[extinction_probability > 1] <- 1

  colonisation_probability <- ColProb(
    patch_area,
    Oij,
    dij,
    alpha,
    ygamma
  )

  colonisation_probability <- .rast(
    assignPatchValues(patch_matrix, colonisation_probability),
    patch
  )

  occupancy <-
    colonisation_probability /
    (
      colonisation_probability +
      extinction_probability -
      colonisation_probability * extinction_probability
    )

  if (all(!is.finite(terra::values(occupancy)))) {
    return(NULL)
  }

  occupancy
}


# Run repeated stochastic realisations of the hybrid model and return their mean.
mean_hybrid_surface <- function(n_runs, ..., verbose = FALSE) {
  runs <- vector("list", n_runs)

  for (i in seq_len(n_runs)) {
    if (verbose) {
      message("Simulation ", i, " of ", n_runs)
    }

    runs[[i]] <- OptimMetaPopHybrid(...)
  }

  runs <- Filter(Negate(is.null), runs)

  if (length(runs) == 0) {
    return(NULL)
  }

  mean(terra::rast(runs), na.rm = TRUE)
}


# Pearson correlation between corresponding cells of two occupancy surfaces.
pearson_surface_cor <- function(true_surface, predicted_surface) {
  values <- cbind(
    terra::values(true_surface),
    terra::values(predicted_surface)
  )

  values <- values[stats::complete.cases(values), , drop = FALSE]

  if (nrow(values) < 3) {
    return(NA_real_)
  }

  stats::cor(values[, 1], values[, 2], method = "pearson")
}


# Warren's I similarity statistic for two non-negative raster surfaces.
warrens_I <- function(surface_1, surface_2) {
  values <- cbind(
    terra::values(surface_1),
    terra::values(surface_2)
  )

  values <- values[stats::complete.cases(values), , drop = FALSE]

  if (nrow(values) == 0) {
    return(NA_real_)
  }

  values[, 1] <- pmax(values[, 1], 0)
  values[, 2] <- pmax(values[, 2], 0)

  total_1 <- sum(values[, 1])
  total_2 <- sum(values[, 2])

  if (total_1 == 0 || total_2 == 0) {
    return(NA_real_)
  }

  prob_1 <- values[, 1] / total_1
  prob_2 <- values[, 2] / total_2

  sum(sqrt(prob_1 * prob_2))
}


# ==============================================================================
# 3. REFERENCE SIMULATION
# ==============================================================================

true_pars <- c(
  p1 = 1.5,
  p2 = 1.0,
  p3 = 2.0,
  p4 = -5.0,
  mu = 0.01,
  x = 1.0,
  alpha = 10.0,
  ygamma = 2.0,
  aos = 0.1
)

# Analysis settings used throughout Box 2.
timesteps <- 500
min_density <- 0.01
max_density <- 20
n_reference_runs <- 50
n_objective_runs <- 10
n_optimizer_runs <- 20
n_validation_runs <- 50

# Prepare habitat predictors.
lcover <- terra::rast(landcover)
vhgt <- terra::rast(veghgt)

broadwood <- lcover * 0
conifwood <- lcover * 0

broadwood[lcover == 1] <- 1
conifwood[lcover == 2] <- 1

msk <- 0.5 * (broadwood + conifwood)
msk[msk == 0] <- NA

# Generate the reference occupancy surface.
averageTrue <- generate_reference_surface(
  habclass1 = conifwood,
  habclass2 = broadwood,
  habcont = vhgt,
  msk = msk,
  p1 = true_pars["p1"],
  p2 = true_pars["p2"],
  p3 = true_pars["p3"],
  p4 = true_pars["p4"],
  mu = true_pars["mu"],
  x = true_pars["x"],
  alpha = true_pars["alpha"],
  ygamma = true_pars["ygamma"],
  aos = true_pars["aos"],
  timesteps = timesteps,
  minden = min_density,
  maxden = max_density,
  asprob = TRUE,
  n_runs = n_reference_runs
)


# ==============================================================================
# 4. BAYESIAN OPTIMISATION
# ==============================================================================

parameter_names <- c(
  "p1", "p2", "p3", "p4",
  "mu", "x", "alpha", "ygamma"
)

full_bounds <- list(
  p1 = c(-3, 3),
  p2 = c(-3, 3),
  p3 = c(-3, 3),
  p4 = c(-8, -1),
  mu = c(0.005, 0.2),
  x = c(0.3, 2),
  alpha = c(5, 15),
  ygamma = c(0.5, 4)
)


# Objective function passed to rBayesianOptimization.
# Warren's I is maximised against the complete reference occupancy surface.
objective_function <- function(
    habclass1 = conifwood,
    habclass2 = broadwood,
    habcont = vhgt,
    msk = msk,
    p1,
    p2,
    p3,
    p4,
    mu,
    x,
    alpha,
    ygamma,
    aos = true_pars["aos"],
    n_runs = n_objective_runs
 ) {
  average_predicted <- mean_hybrid_surface(
    n_runs = n_runs,
    habclass1 = habclass1,
    habclass2 = habclass2,
    habcont = habcont,
    msk = msk,
    p1 = p1,
    p2 = p2,
    p3 = p3,
    p4 = p4,
    mu = mu,
    x = x,
    alpha = alpha,
    ygamma = ygamma,
    aos = aos,
    timesteps = timesteps,
    minden = min_density,
    maxden = max_density
  )

  if (is.null(average_predicted)) {
    return(list(Score = 0.001, Pred = 0))
  }

  similarity <- warrens_I(average_predicted, averageTrue)

  if (!is.finite(similarity)) {
    return(list(Score = 0.001, Pred = 0))
  }

  list(Score = similarity, Pred = 0)
}

# Run the optimiser independently.
opt_results <- vector("list", n_optimizer_runs)

for (i in seq_len(n_optimizer_runs)) {
  message("Starting optimisation run ", i, " of ", n_optimizer_runs)

  opt_res <- rBayesianOptimization::BayesianOptimization(
    objective_function,
    bounds = full_bounds,
    init_points = 50,
    n_iter = 50,
    acq = "ucb",
    kappa = 2.576,
    eps = 0,
    kernel = list(type = "matern", nu = 2.5),
    verbose = TRUE
  )

  # Store the best parameter combination and its objective value.
  opt_results[[i]] <- c(
    opt_res$Best_Par,
    value = opt_res$Best_Value
  )
}

# Convert optimisation results to a single data frame.
opt_df <- purrr::imap_dfr(
  opt_results,
  function(result, run_id) {
    tibble::as_tibble_row(as.list(result)) |>
      dplyr::mutate(run_id = run_id, .before = 1)
  }
)

# Save optimisation results.
utils::write.csv(
  opt_df,
  file.path(results_dir, "bayesian_optimisation_parameter_estimates.csv"),
  row.names = FALSE
)




# ==============================================================================
# 5. SELECT THE HIGHEST-SCORING PARAMETER COMBINATION
# ==============================================================================

best_run <- opt_df[which.max(opt_df$value), , drop = FALSE]

best_pars <- setNames(
  as.numeric(best_run[1, parameter_names]),
  parameter_names
)

message(
  "Highest-scoring optimisation run: ", best_run$run_id,
  " (Warren's I = ", round(best_run$value, 6), ")"
)

print(best_pars)


# ==============================================================================
# 6. VALIDATE THE HIGHEST-SCORING PARAMETERISATION
# ==============================================================================

averageBayesPred <- mean_hybrid_surface(
  n_runs = n_validation_runs,
  habclass1 = conifwood,
  habclass2 = broadwood,
  habcont = vhgt,
  msk = msk,
  p1 = best_pars["p1"],
  p2 = best_pars["p2"],
  p3 = best_pars["p3"],
  p4 = best_pars["p4"],
  mu = best_pars["mu"],
  x = best_pars["x"],
  alpha = best_pars["alpha"],
  ygamma = best_pars["ygamma"],
  aos = true_pars["aos"],
  timesteps = timesteps,
  minden = min_density,
  maxden = max_density,
  verbose = TRUE
)

if (is.null(averageBayesPred)) {
  stop("All validation simulations failed.")
}

surface_comparison <- data.frame(
  pearson_r = pearson_surface_cor(averageTrue, averageBayesPred),
  warrens_I = warrens_I(averageTrue, averageBayesPred)
)

print(surface_comparison)

utils::write.csv(
  surface_comparison,
  file.path(results_dir, "surface_comparison_metrics.csv"),
  row.names = FALSE
)


# ==============================================================================
# 7. FIGURE 1: REFERENCE AND BAYESIAN-PREDICTED OCCUPANCY SURFACES
# ==============================================================================
# Panel a = reference occupancy surface generated using known parameters.
# Panel b = mean occupancy surface from 50 simulations using the highest-scoring
#           Bayesian optimisation parameter combination.
# ==============================================================================

figure1_file <- file.path(
  figures_dir,
  "figure1_reference_vs_bayesian_occupancy.png"
)

occupancy_opacity <- 0.88

if (!terra::compareGeom(averageBayesPred, averageTrue, stopOnError = FALSE)) {
  stop("Predicted and reference surfaces do not have matching raster geometry.")
}

if (is.na(terra::crs(averageTrue)) || terra::crs(averageTrue) == "") {
  stop("The reference occupancy surface must have a valid CRS.")
}

if (terra::is.lonlat(averageTrue)) {
  warning(
    "The map scale bar assumes projected coordinates in metres, but the raster ",
    "appears to use longitude/latitude coordinates."
  )
}

# Drawing copies: zero/NA cells are made transparent through separate alpha rasters.
reference_draw <- averageTrue
bayes_draw <- averageBayesPred

reference_alpha <- terra::ifel(
  !is.na(reference_draw) & reference_draw > 0,
  occupancy_opacity,
  0
)

bayes_alpha <- terra::ifel(
  !is.na(bayes_draw) & bayes_draw > 0,
  occupancy_opacity,
  0
)

reference_draw[is.na(reference_draw)] <- 0
bayes_draw[is.na(bayes_draw)] <- 0

# Map extent.
map_ext <- terra::ext(averageTrue)
xmin <- map_ext$xmin
xmax <- map_ext$xmax
ymin <- map_ext$ymin
ymax <- map_ext$ymax
map_width <- xmax - xmin
map_height <- ymax - ymin

# Download and cache Natural Earth land data locally.
ne_zip <- file.path(cache_dir, "ne_10m_land.zip")
ne_shp <- file.path(cache_dir, "ne_10m_land.shp")
ne_url <- "https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_land.zip"

if (!file.exists(ne_shp)) {
  message("Downloading Natural Earth 1:10m land basemap...")
  utils::download.file(ne_url, ne_zip, mode = "wb", quiet = FALSE)
  utils::unzip(ne_zip, exdir = cache_dir)
}

if (!file.exists(ne_shp)) {
  stop("Natural Earth land shapefile was not extracted successfully.")
}

world_land <- terra::vect(ne_shp)
world_land <- terra::project(world_land, terra::crs(averageTrue))
world_land <- terra::crop(world_land, map_ext)

# Study-area outline.
study_mask <- terra::ifel(!is.na(lcover), 1, NA)
study_outline <- terra::as.polygons(
  study_mask,
  dissolve = TRUE,
  na.rm = TRUE
)

# Figure colours.
sea_colour <- "grey75"
land_colour <- "grey91"
coast_colour <- "grey45"
map_colours <- grDevices::hcl.colors(120, palette = "Viridis")
common_range <- c(0, 1)


draw_basemap <- function() {
  graphics::plot(
    NA,
    xlim = c(xmin, xmax),
    ylim = c(ymin, ymax),
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = "",
    xaxs = "i",
    yaxs = "i"
  )

  graphics::rect(
    xmin,
    ymin,
    xmax,
    ymax,
    col = sea_colour,
    border = NA
  )

  if (nrow(world_land) > 0) {
    terra::plot(
      world_land,
      add = TRUE,
      col = land_colour,
      border = coast_colour,
      lwd = 0.65
    )
  }
}


add_panel_label <- function(label) {
  x <- xmin + 0.025 * map_width
  y <- ymax - 0.025 * map_height

  graphics::rect(
    x - 0.008 * map_width,
    y - 0.070 * map_height,
    x + 0.060 * map_width,
    y + 0.015 * map_height,
    col = grDevices::adjustcolor("white", alpha.f = 0.88),
    border = NA
  )

  graphics::text(
    x,
    y,
    labels = label,
    adj = c(0, 1),
    font = 2,
    cex = 1.4
  )
}


add_north_arrow <- function() {
  x <- xmin + 0.86 * map_width
  y0 <- ymin + 0.12 * map_height
  y1 <- ymin + 0.25 * map_height

  graphics::rect(
    x - 0.038 * map_width,
    y0 - 0.035 * map_height,
    x + 0.038 * map_width,
    y1 + 0.065 * map_height,
    col = grDevices::adjustcolor("white", alpha.f = 0.80),
    border = NA
  )

  graphics::text(
    x,
    y1 + 0.032 * map_height,
    "N",
    font = 2,
    cex = 0.95
  )

  graphics::arrows(
    x0 = x,
    y0 = y0,
    x1 = x,
    y1 = y1,
    length = 0.10,
    lwd = 1.6
  )
}


add_scale_bar <- function() {
  # Assumes map coordinates are expressed in metres.
  total_length <- 5000
  segment_length <- 2500

  x0 <- xmin + 0.12 * map_width
  y0 <- ymin + 0.13 * map_height
  bar_height <- 0.030 * map_height

  graphics::rect(
    x0 - 0.025 * map_width,
    y0 - 0.085 * map_height,
    x0 + total_length + 0.025 * map_width,
    y0 + 0.060 * map_height,
    col = grDevices::adjustcolor("white", alpha.f = 0.83),
    border = NA
  )

  graphics::rect(
    x0,
    y0,
    x0 + segment_length,
    y0 + bar_height,
    col = "grey15",
    border = "grey10",
    lwd = 1
  )

  graphics::rect(
    x0 + segment_length,
    y0,
    x0 + total_length,
    y0 + bar_height,
    col = "white",
    border = "grey10",
    lwd = 1
  )

  graphics::segments(
    c(x0, x0 + segment_length, x0 + total_length),
    y0,
    c(x0, x0 + segment_length, x0 + total_length),
    y0 - 0.018 * map_height,
    lwd = 1.1
  )

  graphics::text(
    c(x0, x0 + segment_length, x0 + total_length),
    y0 - 0.033 * map_height,
    labels = c("0", "2.5", "5 km"),
    cex = 0.82,
    adj = c(0.5, 1)
  )
}


draw_map_panel <- function(raster, alpha_raster, label, cartography = FALSE) {
  graphics::par(
    mar = c(0.05, 0.05, 0.05, 0.05),
    xaxs = "i",
    yaxs = "i"
  )

  draw_basemap()

  terra::plot(
    raster,
    add = TRUE,
    range = common_range,
    fill_range = TRUE,
    col = map_colours,
    alpha = alpha_raster,
    colNA = NA,
    background = NULL,
    axes = FALSE,
    legend = FALSE,
    box = FALSE,
    reset = FALSE
  )

  terra::plot(
    study_outline,
    add = TRUE,
    col = NA,
    border = "grey20",
    lwd = 0.85
  )

  if (cartography) {
    add_north_arrow()
    add_scale_bar()
  }

  add_panel_label(label)
}


draw_colour_bar <- function() {
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::plot.window(
    xlim = c(0, 1),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i"
  )

  left <- 0.32
  right <- 0.68
  bottom <- 0.48
  top <- 0.70

  breaks <- seq(0, 1, length.out = length(map_colours) + 1)

  for (i in seq_along(map_colours)) {
    x1 <- left + (right - left) * breaks[i]
    x2 <- left + (right - left) * breaks[i + 1]
    graphics::rect(x1, bottom, x2, top, col = map_colours[i], border = NA)
  }

  graphics::rect(
    left,
    bottom,
    right,
    top,
    border = "grey20",
    lwd = 0.8
  )

  ticks <- seq(0, 1, by = 0.2)
  tick_x <- left + (right - left) * ticks

  graphics::segments(
    tick_x,
    bottom,
    tick_x,
    bottom - 0.06,
    lwd = 0.7
  )

  graphics::text(
    tick_x,
    bottom - 0.10,
    labels = sprintf("%.1f", ticks),
    cex = 0.72,
    adj = c(0.5, 1)
  )

  graphics::text(
    0.5,
    0.91,
    "Probability of presence",
    cex = 0.80
  )
}


draw_reference_bayes_figure <- function() {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  graphics::layout(
    matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE),
    widths = c(1, 1),
    heights = c(1, 0.11)
  )

  # Panel a: reference occupancy surface.
  draw_map_panel(
    raster = reference_draw,
    alpha_raster = reference_alpha,
    label = "a",
    cartography = TRUE
  )

  # Panel b: Bayesian-optimised occupancy surface.
  draw_map_panel(
    raster = bayes_draw,
    alpha_raster = bayes_alpha,
    label = "b",
    cartography = FALSE
  )

  draw_colour_bar()
}

if (interactive()) {
  grDevices::dev.new(width = 9, height = 4.4)
  draw_reference_bayes_figure()
}

grDevices::png(
  filename = figure1_file,
  width = 3600,
  height = 1700,
  units = "px",
  res = 400,
  bg = "white"
)

draw_reference_bayes_figure()
grDevices::dev.off()
message("Saved Figure 1: ", figure1_file)


# ==============================================================================
# 8. FIGURE 2: PARAMETER ESTIMATES ACROSS OPTIMISATION RUNS
# ==============================================================================

# Use the same generating values and optimisation bounds as the model itself so
# that the figure cannot silently become inconsistent with the analysis.
true_values <- tibble::tibble(
  parameter = parameter_names,
  true_value = as.numeric(true_pars[parameter_names])
)

parameter_bounds <- tibble::tibble(
  parameter = parameter_names,
  lower = purrr::map_dbl(full_bounds[parameter_names], 1),
  upper = purrr::map_dbl(full_bounds[parameter_names], 2)
)

# Convert the optimisation results directly from opt_df.
plot_data <- opt_df |>
  dplyr::select(run_id, dplyr::all_of(parameter_names)) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(parameter_names),
    names_to = "parameter",
    values_to = "estimate"
  ) |>
  dplyr::mutate(
    parameter = factor(parameter, levels = parameter_names)
  )

true_values <- true_values |>
  dplyr::mutate(
    parameter = factor(parameter, levels = parameter_names)
  )

bounds_long <- parameter_bounds |>
  tidyr::pivot_longer(
    cols = c(lower, upper),
    names_to = "bound",
    values_to = "bound_value"
  ) |>
  dplyr::mutate(
    parameter = factor(parameter, levels = parameter_names)
  )

# Wrapped facet labels with model symbols.
parameter_labels <- c(
  p1 = "Response to conifer<br>woodland<br><i>p</i><sub>1</sub>",
  p2 = "Response to broadleaf<br>woodland<br><i>p</i><sub>2</sub>",
  p3 = "Response to vegetation<br>height<br><i>p</i><sub>3</sub>",
  p4 = "SDM intercept<br><i>p</i><sub>4</sub>",
  mu = "Baseline extinction<br>&mu;",
  x = "Extinction scaling<br><i>x</i>",
  alpha = "Dispersal<br>&alpha;",
  ygamma = "Colonisation scaling<br>&gamma;"
)

# Reproducible jitter positions.
set.seed(123)

fig2 <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = 1, y = estimate)
) +
  # Bounds imposed during Bayesian optimisation.
  ggplot2::geom_hline(
    data = bounds_long,
    ggplot2::aes(yintercept = bound_value),
    inherit.aes = FALSE,
    colour = "grey60",
    linetype = "dotted",
    linewidth = 0.45
  ) +
  # Distribution across independent optimisation runs.
  ggplot2::geom_boxplot(
    width = 0.30,
    fill = "white",
    outlier.shape = NA,
    linewidth = 0.55
  ) +
  # Individual estimates from all 20 optimisation runs.
  ggplot2::geom_jitter(
    width = 0.10,
    height = 0,
    size = 2.0,
    alpha = 0.75
  ) +
  # Known generating value for each parameter.
  ggplot2::geom_hline(
    data = true_values,
    ggplot2::aes(yintercept = true_value),
    inherit.aes = FALSE,
    colour = "#D55E00",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  ggplot2::facet_wrap(
    ~parameter,
    scales = "free_y",
    ncol = 4,
    labeller = ggplot2::as_labeller(parameter_labels)
  ) +
  ggplot2::scale_x_continuous(
    breaks = NULL,
    name = NULL
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.05, 0.08))
  ) +
  ggplot2::labs(
    y = "Estimated parameter value"
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggtext::element_markdown(
      size = 8.5,
      lineheight = 1.05,
      margin = ggplot2::margin(t = 3, b = 4)
    ),
    axis.title.y = ggplot2::element_text(
      size = 10,
      margin = ggplot2::margin(r = 7)
    ),
    axis.text.y = ggplot2::element_text(size = 8.5),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(1.2, "lines"),
    panel.spacing.y = grid::unit(1.2, "lines"),
    plot.margin = ggplot2::margin(t = 6, r = 8, b = 6, l = 6)
  )

if (interactive()) {
  print(fig2)
}

# Vector version.
ggplot2::ggsave(
  filename = file.path(figures_dir, "figure2_parameter_estimates.pdf"),
  plot = fig2,
  width = 180,
  height = 110,
  units = "mm",
  device = grDevices::cairo_pdf
)

# High-resolution raster version.
ggplot2::ggsave(
  filename = file.path(figures_dir, "figure2_parameter_estimates.tiff"),
  plot = fig2,
  width = 180,
  height = 110,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

message("Saved Figure 2 to: ", figures_dir)


# ==============================================================================
# 9. REPRODUCIBILITY INFORMATION
# ==============================================================================

writeLines(
  capture.output(utils::sessionInfo()),
  file.path(results_dir, "session_info.txt")
)

message("Analysis complete.")
