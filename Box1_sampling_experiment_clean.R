# ==============================================================================
# BOX 1: SAMPLING AND RECONSTRUCTION IN A HYBRID SDM-METAPOPULATION MODEL
# ==============================================================================
# Purpose
#   Examine how sampling intensity and minimum spatial separation affect the
#   ability of fitmetapopsdm() to reconstruct a known metapopulation occupancy
#   surface and recover fitted metapopulation parameters.
#
# Workflow
#   1. Construct habitat suitability from the PopScape example landscape.
#   2. Simulate and average metapopulation occupancy to create a reference surface.
#   3. Draw repeated presence-absence samples across combinations of sample size
#      and minimum spatial separation.
#   4. Fit the hybrid SDM-metapopulation model to each sample.
#   5. Summarise numerical convergence, parameter recovery and reconstruction
#      accuracy using Warren's I and Pearson's r.
#   6. Save analysis outputs and generate Figure 1.
#
# Notes
#   - Sampled occupancy probabilities are converted to presence-absence at 0.5.
#   - Surface metrics compare the reference surface with the mean reconstructed
#     surface across numerically converged fits within each sampling configuration.
#   - A fit is classified as converged when the fitted parameters used in the
#     convergence check are finite; this is a numerical, not biological, criterion.
#   - alpha and aos are fixed at their generating values during fitting, so
#     parameter-recovery summaries focus on mu, x and ygamma.
# ==============================================================================

library(terra)
library(PopScape)

# ==============================================================================
# 1. SETTINGS
# ==============================================================================

set.seed(123)

n_true_runs <- 50
n_fit_reps  <- 50

sample_sizes  <- c(50, 100, 300, 400)
min_distances <- c(0, 50, 100, 200)  # metres for this projected landscape

# At least one observation from each valid land-cover class is required so the
# categorical predictor can be fitted. Remaining points are sampled randomly.
min_per_class <- 1
equal_classes <- FALSE

# Figure settings
make_figure        <- TRUE
show_sample_points <- TRUE
preview_figure     <- FALSE
figure_file        <- "Figure_1_occupancy_reconstruction.png"

# Project-relative output directories.
output_dir <- file.path("outputs", "box1")
cache_dir  <- file.path("cache", "natural_earth")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. LANDSCAPE AND HABITAT SUITABILITY
# ==============================================================================

lcover <- rast(landcover)
vhgt   <- rast(veghgt)

if (!compareGeom(lcover, vhgt, stopOnError = FALSE)) {
  stop("landcover and veghgt do not have matching raster geometry.")
}
if (is.lonlat(lcover)) {
  stop("This experiment uses minDist in metres; project the landscape before running.")
}

broadwood <- lcover * 0
broadwood[lcover == 1] <- 1

conifwood <- lcover * 0
conifwood[lcover == 2] <- 1

habitat_mask <- 0.5 * (broadwood + conifwood)
habitat_mask[habitat_mask == 0] <- NA

linear_habitat_suitability <-
  1.5 * conifwood +
  1.0 * broadwood +
  2.0 * log(vhgt) -
  5

habitat_suitability <- mask(1 / (1 + exp(-linear_habitat_suitability)), habitat_mask)

# ==============================================================================
# 3. METAPOPULATION SIMULATION
# ==============================================================================

metapop_params <- list(
  mu        = 0.01,
  x         = 1,
  alpha     = 10,
  ygamma    = 2,
  aos       = 0.1,
  timesteps = 500,
  minden    = 0.01,
  maxden    = 20,
  asprob    = TRUE
)

true_occupancy_runs <- vector("list", n_true_runs)

for (i in seq_len(n_true_runs)) {
  message("Running MetaPopSim ", i, " / ", n_true_runs)
  
  true_occupancy_runs[[i]] <- MetaPopSim(
    habitat_suitability,
    mu        = metapop_params$mu,
    x         = metapop_params$x,
    alpha     = metapop_params$alpha,
    ygamma    = metapop_params$ygamma,
    aos       = metapop_params$aos,
    timesteps = metapop_params$timesteps,
    minden    = metapop_params$minden,
    maxden    = metapop_params$maxden,
    asprob    = metapop_params$asprob
  )
}

# The simulations are averaged first. Every sampling configuration below uses
# this same probability-of-occupancy surface.
average_true_occupancy <- mean(rast(true_occupancy_runs))

# ==============================================================================
# 4. SPATIALLY CONSTRAINED SUBSAMPLING
# ==============================================================================
# Based closely on PopScape::subsample(), with three additions:
#   minDist        minimum Euclidean distance between sampled locations;
#   min_per_class  minimum representation of each valid land-cover class;
#   equal_classes  optional approximately balanced sampling across classes.
#
# As in the PopScape sampling logic used for this experiment, cells with exactly
# zero occupancy are unavailable for sampling and sampled probability values are
# converted to presence/absence using a threshold of 0.5.
# ==============================================================================

subsample_spatial <- function(occo,
                              landcover,
                              n,
                              minDist = 0,
                              min_per_class = 1,
                              equal_classes = FALSE) {
  if (n <= 0)
    stop("n must be greater than zero.")
  if (minDist < 0)
    stop("minDist must be >= 0.")
  if (min_per_class < 0)
    stop("min_per_class must be >= 0.")
  
  # Do not modify the supplied raster outside this function.
  occo <- occo * 1
  occo[occo == 0] <- NA
  
  landcovers <- mask(landcover, occo)
  suithabs <- unique(as.vector(landcovers))
  suithabs <- suithabs[!is.na(suithabs)]
  n_classes <- length(suithabs)
  
  if (n_classes < 1)
    stop("No valid land-cover classes found in the occupancy area.")
  if (min_per_class * n_classes > n) {
    stop("Sample size is too small to satisfy min_per_class for all land-cover classes.")
  }
  
  # Reproduce the PopScape candidate-density calculation within suitable classes.
  hab <- landcover * NA
  for (hab_class in suithabs)
    hab[landcover == hab_class] <- 1
  landcover_sample_area <- mask(landcover, hab)
  
  v <- as.vector(landcover_sample_area)
  mu_candidate <- length(v) / sum(!is.na(v))
  n_candidates <- trunc(n * mu_candidate * 2)
  if (minDist > 0)
    n_candidates <- max(n_candidates, n * 20)
  
  e <- ext(occo)
  candidates <- data.frame(
    x = runif(n_candidates, e$xmin, e$xmax),
    y = runif(n_candidates, e$ymin, e$ymax)
  )
  
  candidates$pa <- terra::extract(occo, candidates)[, 2]
  candidates <- candidates[!is.na(candidates$pa), , drop = FALSE]
  candidates$pa <- ifelse(candidates$pa < 0.5, 0, 1)
  
  candidates$landcover <- terra::extract(landcover, candidates[, c("x", "y"), drop = FALSE])[, 2]
  candidates <- candidates[!is.na(candidates$landcover), , drop = FALSE]
  
  # Add extra candidates from rare classes where the general random pool is thin.
  add_class_candidates <- function(hab_class, n_required) {
    class_raster <- landcovers
    class_raster[class_raster != hab_class] <- NA
    
    class_points <- terra::spatSample(
      class_raster,
      size = max(100, n_required * 20),
      method = "random",
      na.rm = TRUE,
      xy = TRUE
    )
    
    if (is.null(class_points) ||
        nrow(class_points) == 0)
      return(NULL)
    
    class_points <- class_points[, c("x", "y"), drop = FALSE]
    class_points$pa <- terra::extract(occo, class_points)[, 2]
    class_points <- class_points[!is.na(class_points$pa), , drop = FALSE]
    class_points$pa <- ifelse(class_points$pa < 0.5, 0, 1)
    class_points$landcover <- hab_class
    class_points
  }
  
  if (min_per_class > 0 || equal_classes) {
    for (hab_class in suithabs) {
      n_existing <- sum(candidates$landcover == hab_class)
      target_pool <- max(min_per_class, 10)
      
      if (n_existing < target_pool) {
        extra <- add_class_candidates(hab_class, target_pool)
        if (!is.null(extra))
          candidates <- rbind(candidates, extra)
      }
    }
  }
  
  if (nrow(candidates) == 0)
    stop("No candidate sampling locations were available.")
  candidates <- candidates[sample(seq_len(nrow(candidates))), , drop = FALSE]
  
  can_add_point <- function(candidate_index, selected_indices) {
    if (length(selected_indices) == 0 || minDist == 0)
      return(TRUE)
    
    selected_coords <- candidates[selected_indices, c("x", "y"), drop = FALSE]
    d <- sqrt((candidates$x[candidate_index] - selected_coords$x)^2 +
                (candidates$y[candidate_index] - selected_coords$y)^2)
    all(d >= minDist)
  }
  
  # sample(x, length(x)) avoids R's special behaviour when x is a single number.
  shuffle <- function(x) {
    if (length(x) <= 1)
      return(x)
    sample(x, length(x))
  }
  
  select_from_class <- function(hab_class,
                                n_required,
                                selected_indices) {
    if (n_required <= 0)
      return(selected_indices)
    
    class_candidates <- which(candidates$landcover == hab_class)
    class_candidates <- setdiff(class_candidates, selected_indices)
    class_candidates <- shuffle(class_candidates)
    n_added <- 0
    
    for (candidate in class_candidates) {
      if (can_add_point(candidate, selected_indices)) {
        selected_indices <- c(selected_indices, candidate)
        n_added <- n_added + 1
      }
      if (n_added >= n_required)
        break
    }
    selected_indices
  }
  
  # First guarantee the minimum representation of each valid class.
  selected <- integer(0)
  if (min_per_class > 0) {
    for (hab_class in suithabs) {
      n_before <- length(selected)
      selected <- select_from_class(hab_class, min_per_class, selected)
      
      if ((length(selected) - n_before) < min_per_class) {
        stop(
          paste0(
            "Could not obtain ",
            min_per_class,
            " sample(s) from land-cover class ",
            hab_class,
            " while respecting minDist = ",
            minDist,
            "."
          )
        )
      }
    }
  }
  
  # Optional balanced design. This is not used in the manuscript experiment
  # because equal_classes = FALSE, but is retained as a sampling option.
  if (equal_classes) {
    targets <- rep(floor(n / n_classes), n_classes)
    remainder <- n - sum(targets)
    if (remainder > 0)
      targets[seq_len(remainder)] <- targets[seq_len(remainder)] + 1
    if (any(targets < min_per_class)) {
      stop("equal_classes = TRUE is incompatible with min_per_class at this sample size.")
    }
    names(targets) <- suithabs
    
    for (j in seq_along(suithabs)) {
      hab_class <- suithabs[j]
      current <- sum(candidates$landcover[selected] == hab_class)
      selected <- select_from_class(hab_class, targets[j] - current, selected)
      
      current_after <- sum(candidates$landcover[selected] == hab_class)
      if (current_after < targets[j]) {
        stop(
          paste0(
            "Could not achieve the equal-class target for land-cover class ",
            hab_class,
            " with minDist = ",
            minDist,
            "."
          )
        )
      }
    }
  }
  
  # Fill the remainder randomly while respecting minDist.
  remaining_candidates <- setdiff(seq_len(nrow(candidates)), selected)
  remaining_candidates <- shuffle(remaining_candidates)
  
  for (candidate in remaining_candidates) {
    if (length(selected) >= n)
      break
    if (can_add_point(candidate, selected))
      selected <- c(selected, candidate)
  }
  
  if (length(selected) < n) {
    stop(
      paste0(
        "Only ",
        length(selected),
        " points could be obtained; ",
        n,
        " requested with minDist = ",
        minDist,
        "."
      )
    )
  }
  
  dfout <- candidates[selected[seq_len(n)], c("x", "y", "pa"), drop = FALSE]
  
  if (length(unique(dfout$pa)) < 2) {
    stop("Final sample contains only one presence/absence class.")
  }
  
  final_lc <- terra::extract(landcover, dfout[, c("x", "y"), drop = FALSE])[, 2]
  final_counts <- table(factor(final_lc, levels = suithabs))
  if (any(final_counts < min_per_class)) {
    stop("Final sample does not satisfy min_per_class.")
  }
  
  dfout
}

# ==============================================================================
# 5. SURFACE-ACCURACY FUNCTIONS
# ==============================================================================

warrens_I <- function(surface_1, surface_2) {
  vals <- cbind(terra::values(surface_1), terra::values(surface_2))
  vals <- vals[complete.cases(vals), , drop = FALSE]
  if (nrow(vals) == 0)
    return(NA_real_)
  
  vals[, 1] <- pmax(vals[, 1], 0)
  vals[, 2] <- pmax(vals[, 2], 0)
  total_1 <- sum(vals[, 1])
  total_2 <- sum(vals[, 2])
  if (total_1 == 0 || total_2 == 0)
    return(NA_real_)
  
  p1 <- vals[, 1] / total_1
  p2 <- vals[, 2] / total_2
  sum(sqrt(p1 * p2))
}

pearson_surface_cor <- function(true_surface, predicted_surface) {
  vals <- cbind(terra::values(true_surface),
                terra::values(predicted_surface))
  vals <- vals[complete.cases(vals), , drop = FALSE]
  if (nrow(vals) < 3)
    return(NA_real_)
  if (sd(vals[, 1]) == 0 || sd(vals[, 2]) == 0)
    return(NA_real_)
  stats::cor(vals[, 1], vals[, 2], method = "pearson")
}

# ==============================================================================
# 6. SAMPLING CONFIGURATIONS AND STORAGE
# ==============================================================================

sampling_configs <- expand.grid(
  n_samples = sample_sizes,
  minDist = min_distances,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
sampling_configs$config_id <- seq_len(nrow(sampling_configs))
sampling_configs$min_per_class <- min_per_class
sampling_configs$equal_classes <- equal_classes

subsamples <- vector("list", nrow(sampling_configs))
predicted_surfaces <- vector("list", nrow(sampling_configs))
predicted_surfaces_converged <- vector("list", nrow(sampling_configs))
parameter_results <- list()
surface_results <- vector("list", nrow(sampling_configs))

names(subsamples) <- names(predicted_surfaces) <- names(predicted_surfaces_converged) <-
  paste0("config_", sampling_configs$config_id)

# These parameters must be finite for a fit to be called numerically converged.
# alpha and aos are fixed in this experiment but are included defensively because
# they are returned in fitted$params.
convergence_parameters <- c("mu", "x", "alpha", "ygamma", "aos")

# Helper used because failed rows and fitted-parameter rows have different columns.
bind_fill <- function(x) {
  all_names <- unique(unlist(lapply(x, names)))
  x <- lapply(x, function(df) {
    missing_names <- setdiff(all_names, names(df))
    for (nm in missing_names)
      df[[nm]] <- NA
    df[, all_names, drop = FALSE]
  })
  out <- do.call(rbind, x)
  row.names(out) <- NULL
  out
}

# ==============================================================================
# 7. SAMPLING AND MODEL FITTING
# ==============================================================================
# Failures are experimental outcomes and are recorded rather than stopping the
# full experiment. Sampling/model errors are captured on the original attempt;
# failed stochastic operations are not rerun simply to retrieve an error message.
# ==============================================================================

for (config in seq_len(nrow(sampling_configs))) {
  n_samples <- sampling_configs$n_samples[config]
  minDist   <- sampling_configs$minDist[config]
  
  message(
    "\nConfiguration ",
    config,
    " / ",
    nrow(sampling_configs),
    ": n = ",
    n_samples,
    ", minDist = ",
    minDist,
    " m"
  )
  
  config_subsamples <- vector("list", n_fit_reps)
  config_predictions <- vector("list", n_fit_reps)
  config_converged_predictions <- vector("list", n_fit_reps)
  
  for (rep in seq_len(n_fit_reps)) {
    message("  Replicate ", rep, " / ", n_fit_reps)
    
    metadata <- data.frame(
      config_id = config,
      replicate = rep,
      n_samples = n_samples,
      minDist = minDist,
      min_per_class = min_per_class,
      equal_classes = equal_classes,
      status = NA_character_,
      failure_stage = NA_character_,
      failure_message = NA_character_,
      converged = FALSE,
      stringsAsFactors = FALSE
    )
    
    # --------------------------------------------------------------------------
    # Draw occurrence sample
    # --------------------------------------------------------------------------
    sample_result <- tryCatch(
      list(
        ok = TRUE,
        value = subsample_spatial(
          occo = average_true_occupancy,
          landcover = lcover,
          n = n_samples,
          minDist = minDist,
          min_per_class = min_per_class,
          equal_classes = equal_classes
        ),
        error = NA_character_
      ),
      error = function(e)
        list(
          ok = FALSE,
          value = NULL,
          error = conditionMessage(e)
        )
    )
    
    if (!sample_result$ok) {
      metadata$status <- "failed"
      metadata$failure_stage <- "sampling"
      metadata$failure_message <- sample_result$error
      parameter_results[[length(parameter_results) + 1]] <- metadata
      message("    SAMPLING FAILED: ", sample_result$error)
      next
    }
    
    occdata <- sample_result$value
    config_subsamples[[rep]] <- occdata
    
    # --------------------------------------------------------------------------
    # Fit hybrid SDM-metapopulation model
    # --------------------------------------------------------------------------
    fit_result <- tryCatch(
      list(
        ok = TRUE,
        value = fitmetapopsdm(
          occdata = occdata,
          landcover = lcover,
          conpredictors = log(vhgt),
          alpha = metapop_params$alpha,
          aos = metapop_params$aos,
          bwgt = 0.5,
          tol = 0.001,
          maxiter = 100,
          minden = metapop_params$minden,
          maxden = metapop_params$maxden,
          Bayesian = FALSE
        ),
        error = NA_character_
      ),
      error = function(e)
        list(
          ok = FALSE,
          value = NULL,
          error = conditionMessage(e)
        )
    )
    
    if (!fit_result$ok) {
      metadata$status <- "failed"
      metadata$failure_stage <- "model"
      metadata$failure_message <- fit_result$error
      parameter_results[[length(parameter_results) + 1]] <- metadata
      message("    MODEL FIT FAILED: ", fit_result$error)
      next
    }
    
    fitted <- fit_result$value
    fitted_params <- as.data.frame(as.list(fitted$params), check.names = FALSE)
    
    # Store a returned prediction even if the parameterisation later fails the
    # numerical-convergence check; it is retained only for audit/reference.
    predicted <- fitted$pocco
    if (is.null(predicted)) {
      metadata$status <- "failed"
      metadata$failure_stage <- "model"
      metadata$failure_message <- "fitmetapopsdm() returned no pocco surface."
      parameter_results[[length(parameter_results) + 1]] <- cbind(metadata, fitted_params)
      next
    }
    
    predicted[is.na(predicted)] <- 0
    predicted <- mask(predicted, lcover)
    config_predictions[[rep]] <- predicted
    
    # --------------------------------------------------------------------------
    # Numerical convergence classification
    # --------------------------------------------------------------------------
    missing_params <- setdiff(convergence_parameters, names(fitted_params))
    finite_params <- FALSE
    
    if (length(missing_params) == 0) {
      param_values <- as.numeric(unlist(fitted_params[1, convergence_parameters, drop = FALSE], use.names = FALSE))
      finite_params <- all(is.finite(param_values))
    }
    
    if (finite_params) {
      metadata$status <- "success"
      metadata$converged <- TRUE
      config_converged_predictions[[rep]] <- predicted
    } else {
      metadata$status <- "failed"
      metadata$failure_stage <- "non_convergence"
      metadata$failure_message <- if (length(missing_params) > 0) {
        paste("Missing fitted parameter(s):",
              paste(missing_params, collapse = ", "))
      } else {
        "One or more fitted parameters were non-finite (NA, NaN, Inf or -Inf)."
      }
      message("    NON-CONVERGED: ", metadata$failure_message)
    }
    
    parameter_results[[length(parameter_results) + 1]] <- cbind(metadata, fitted_params)
  }
  
  subsamples[[config]] <- config_subsamples
  
  # Retain all returned predictions for audit; manuscript summaries use converged fits only.
  all_predictions <- config_predictions[!vapply(config_predictions, is.null, logical(1))]
  average_all <- if (length(all_predictions) > 0)
    mean(rast(all_predictions))
  else
    NULL
  predicted_surfaces[[config]] <- list(individual = config_predictions, average = average_all)
  
  # Manuscript analyses use only numerically converged fits.
  converged_predictions <- config_converged_predictions[!vapply(config_converged_predictions, is.null, logical(1))]
  
  if (length(converged_predictions) == 0) {
    average_converged <- NULL
    I <- NA_real_
    r <- NA_real_
  } else {
    average_converged <- mean(rast(converged_predictions))
    I <- warrens_I(average_true_occupancy, average_converged)
    r <- pearson_surface_cor(average_true_occupancy, average_converged)
  }
  
  predicted_surfaces_converged[[config]] <- list(individual = config_converged_predictions, average = average_converged)
  
  surface_results[[config]] <- data.frame(
    config_id = config,
    n_samples = n_samples,
    minDist = minDist,
    n_converged = length(converged_predictions),
    warrens_I = I,
    pearson_r = r,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# 8. COMBINE RUN-LEVEL AND CONFIGURATION-LEVEL RESULTS
# ==============================================================================

parameter_df <- bind_fill(parameter_results)
surface_df <- do.call(rbind, surface_results)
row.names(surface_df) <- NULL

warren_df <- surface_df[, c("config_id",
                            "n_samples",
                            "minDist",
                            "n_converged",
                            "warrens_I")]
pearson_df <- surface_df[, c("config_id",
                             "n_samples",
                             "minDist",
                             "n_converged",
                             "pearson_r")]

# Warren's I and Pearson's r are configuration-level metrics, so each value is
# repeated across replicate rows belonging to the same configuration.
parameter_df <- merge(
  parameter_df,
  surface_df[, c("config_id", "warrens_I", "pearson_r")],
  by = "config_id",
  all.x = TRUE,
  sort = FALSE
)
parameter_df <- parameter_df[order(parameter_df$config_id, parameter_df$replicate), , drop = FALSE]

id_columns <- c(
  "config_id",
  "replicate",
  "n_samples",
  "minDist",
  "min_per_class",
  "equal_classes",
  "status",
  "failure_stage",
  "failure_message",
  "warrens_I",
  "pearson_r",
  "converged"
)
parameter_df <- parameter_df[, c(id_columns, setdiff(names(parameter_df), id_columns)), drop = FALSE]
converged_df <- parameter_df[parameter_df$converged, , drop = FALSE]

# ==============================================================================
# 9. CONVERGENCE SUMMARIES
# ==============================================================================

summarise_convergence <- function(df, by) {
  groups <- do.call(interaction, c(df[, by, drop = FALSE], list(drop = TRUE, lex.order = TRUE)))
  
  out <- lapply(split(df, groups), function(x) {
    sampling_failures <- sum(x$failure_stage == "sampling", na.rm = TRUE)
    non_converged <- sum(x$failure_stage == "non_convergence", na.rm = TRUE)
    model_failures <- sum(x$failure_stage == "model", na.rm = TRUE)
    model_attempts <- nrow(x) - sampling_failures
    n_converged <- sum(x$converged, na.rm = TRUE)
    
    cbind(
      x[1, by, drop = FALSE],
      data.frame(
        total_runs = nrow(x),
        model_attempts = model_attempts,
        converged = n_converged,
        non_converged = non_converged,
        model_failures = model_failures,
        sampling_failures = sampling_failures,
        convergence_rate = if (model_attempts > 0)
          n_converged / model_attempts
        else
          NA_real_,
        overall_failure_rate = (non_converged + model_failures + sampling_failures) / nrow(x)
      )
    )
  })
  
  out <- do.call(rbind, out)
  row.names(out) <- NULL
  out
}

convergence_by_n <- summarise_convergence(parameter_df, "n_samples")
convergence_by_distance <- summarise_convergence(parameter_df, "minDist")
convergence_by_config <- summarise_convergence(parameter_df, c("config_id", "n_samples", "minDist"))

# ==============================================================================
# 10. SURFACE-ACCURACY SUMMARY BY SAMPLE SIZE
# ==============================================================================
# Descriptive summaries across the available minimum-distance configurations.
# ==============================================================================

surface_accuracy_by_n <- do.call(rbind, lapply(sample_sizes, function(n_value) {
  x <- surface_df[surface_df$n_samples == n_value, , drop = FALSE]
  w <- x$warrens_I[is.finite(x$warrens_I)]
  p <- x$pearson_r[is.finite(x$pearson_r)]
  
  data.frame(
    n_samples = n_value,
    n_configurations = max(length(w), length(p)),
    mean_warrens_I = if (length(w))
      mean(w)
    else
      NA_real_,
    min_warrens_I = if (length(w))
      min(w)
    else
      NA_real_,
    max_warrens_I = if (length(w))
      max(w)
    else
      NA_real_,
    mean_pearson_r = if (length(p))
      mean(p)
    else
      NA_real_,
    sd_pearson_r = if (length(p) > 1)
      sd(p)
    else
      NA_real_,
    min_pearson_r = if (length(p))
      min(p)
    else
      NA_real_,
    max_pearson_r = if (length(p))
      max(p)
    else
      NA_real_
  )
}))

# ==============================================================================
# 11. PARAMETER-RECOVERY SUMMARY USED IN THE MANUSCRIPT
# ==============================================================================
# alpha and aos are fixed at their known values during fitting, so the recovery
# summary is restricted to mu, x and ygamma. Only numerically converged fits with
# n >= 300 are summarised here. This describes stability/accuracy, not biological
# plausibility; finite but implausible estimates are intentionally retained.
# ==============================================================================

recovery_parameters <- c("mu", "x", "ygamma")
true_parameter_values <- c(mu = metapop_params$mu,
                           x = metapop_params$x,
                           ygamma = metapop_params$ygamma)

high_effort_fits <- converged_df[converged_df$n_samples >= 300, , drop = FALSE]

parameter_recovery_n300plus <- do.call(rbind, lapply(recovery_parameters, function(p) {
  values <- suppressWarnings(as.numeric(high_effort_fits[[p]]))
  values <- values[is.finite(values)]
  q <- if (length(values))
    quantile(values, c(0.25, 0.5, 0.75), names = FALSE)
  else
    rep(NA_real_, 3)
  
  data.frame(
    parameter = p,
    true_value = true_parameter_values[p],
    n_converged = length(values),
    median = q[2],
    q25 = q[1],
    q75 = q[3],
    min = if (length(values))
      min(values)
    else
      NA_real_,
    max = if (length(values))
      max(values)
    else
      NA_real_,
    row.names = NULL
  )
}))

# ==============================================================================
# 12. REPORT KEY RESULTS TO THE CONSOLE
# ==============================================================================

cat("\n--- Convergence by sample size ---\n")
print(convergence_by_n)

cat("\n--- Surface accuracy by sample size ---\n")
print(surface_accuracy_by_n)

cat("\n--- Parameter recovery for converged fits with n >= 300 ---\n")
print(parameter_recovery_n300plus)

# ==============================================================================
# 13. SAVE ANALYSIS OUTPUTS
# ==============================================================================

writeRaster(
  average_true_occupancy,
  file.path(output_dir, "average_true_occupancy.tif"),
  overwrite = TRUE
)

saveRDS(true_occupancy_runs,
        file.path(output_dir, "popscape_true_occupancy_runs.rds"))
saveRDS(subsamples, file.path(output_dir, "popscape_subsamples.rds"))
saveRDS(predicted_surfaces,
        file.path(output_dir, "popscape_predicted_surfaces_all.rds"))
saveRDS(
  predicted_surfaces_converged,
  file.path(output_dir, "popscape_predicted_surfaces_converged.rds")
)

write.csv(parameter_df,
          file.path(output_dir, "popscape_fitted_parameters.csv"),
          row.names = FALSE)
write.csv(
  sampling_configs,
  file.path(output_dir, "popscape_sampling_configurations.csv"),
  row.names = FALSE
)
write.csv(surface_df,
          file.path(output_dir, "popscape_surface_metrics.csv"),
          row.names = FALSE)
write.csv(warren_df,
          file.path(output_dir, "popscape_warrens_I.csv"),
          row.names = FALSE)
write.csv(
  pearson_df,
  file.path(output_dir, "popscape_pearson_surface_correlation.csv"),
  row.names = FALSE
)
write.csv(
  convergence_by_n,
  file.path(output_dir, "popscape_convergence_by_sample_size.csv"),
  row.names = FALSE
)
write.csv(
  convergence_by_distance,
  file.path(output_dir, "popscape_convergence_by_distance.csv"),
  row.names = FALSE
)
write.csv(
  convergence_by_config,
  file.path(output_dir, "popscape_convergence_by_configuration.csv"),
  row.names = FALSE
)
write.csv(
  surface_accuracy_by_n,
  file.path(output_dir, "popscape_surface_accuracy_by_sample_size.csv"),
  row.names = FALSE
)
write.csv(
  parameter_recovery_n300plus,
  file.path(output_dir, "popscape_parameter_recovery_n300plus.csv"),
  row.names = FALSE
)

# ==============================================================================
# 14. FIGURE 1: TRUE AND RECONSTRUCTED OCCUPANCY
# ==============================================================================
# Panel a: true occupancy surface
# Panel b: mean reconstruction, n = 50, minimum distance = 0 m
# Panel c: mean reconstruction, n = 300, minimum distance = 200 m
#
# Predicted surfaces are means across converged fits. If sample points are shown,
# they come from the first converged replicate and illustrate the sampling design;
# they should not be interpreted as the unique sample that generated the mean map.
#
# Natural Earth land polygons are downloaded once and cached locally. The
# occupancy raster uses a cell-level alpha layer so zero/NA cells are transparent
# and the land/sea basemap remains visible.
# ==============================================================================

if (make_figure) {
  config_50_0 <- sampling_configs$config_id[sampling_configs$n_samples == 50 &
                                              sampling_configs$minDist == 0]
  config_300_200 <- sampling_configs$config_id[sampling_configs$n_samples == 300 &
                                                 sampling_configs$minDist == 200]
  
  if (length(config_50_0) != 1 || length(config_300_200) != 1) {
    stop("Could not uniquely identify the two Figure 1 sampling configurations.")
  }
  
  surface_50_0 <- predicted_surfaces_converged[[config_50_0]]$average
  surface_300_200 <- predicted_surfaces_converged[[config_300_200]]$average
  if (is.null(surface_50_0) || is.null(surface_300_200)) {
    stop("One or both Figure 1 configurations have no converged mean prediction.")
  }
  
  occupancy_opacity <- 0.88
  true_draw <- average_true_occupancy
  pred_50_draw <- surface_50_0
  pred_300_draw <- surface_300_200
  
  true_alpha <- ifel(!is.na(true_draw) &
                       true_draw > 0, occupancy_opacity, 0)
  pred_50_alpha <- ifel(!is.na(pred_50_draw) &
                          pred_50_draw > 0, occupancy_opacity, 0)
  pred_300_alpha <- ifel(!is.na(pred_300_draw) &
                           pred_300_draw > 0,
                         occupancy_opacity,
                         0)
  
  true_draw[is.na(true_draw)] <- 0
  pred_50_draw[is.na(pred_50_draw)] <- 0
  pred_300_draw[is.na(pred_300_draw)] <- 0
  
  map_ext <- ext(average_true_occupancy)
  xmin <- map_ext$xmin
  xmax <- map_ext$xmax
  ymin <- map_ext$ymin
  ymax <- map_ext$ymax
  map_width <- xmax - xmin
  map_height <- ymax - ymin
  
  if (is.na(crs(average_true_occupancy)) ||
      crs(average_true_occupancy) == "") {
    stop("average_true_occupancy must have a valid CRS to draw Figure 1.")
  }
  
  # Natural Earth 1:10m land polygons, cached after the first download.
  ne_dir <- cache_dir
  dir.create(ne_dir, showWarnings = FALSE, recursive = TRUE)
  ne_zip <- file.path(ne_dir, "ne_10m_land.zip")
  ne_shp <- file.path(ne_dir, "ne_10m_land.shp")
  ne_url <- "https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_land.zip"
  
  if (!file.exists(ne_shp)) {
    message("Downloading Natural Earth 1:10m land basemap...")
    utils::download.file(ne_url, ne_zip, mode = "wb", quiet = FALSE)
    utils::unzip(ne_zip, exdir = ne_dir)
  }
  if (!file.exists(ne_shp))
    stop("Natural Earth land shapefile was not extracted.")
  
  world_land <- vect(ne_shp)
  world_land <- project(world_land, crs(average_true_occupancy))
  world_land <- crop(world_land, map_ext)
  
  study_mask <- ifel(!is.na(lcover), 1, NA)
  study_outline <- as.polygons(study_mask, dissolve = TRUE, na.rm = TRUE)
  
  sea_colour   <- "grey75"
  land_colour  <- "grey91"
  coast_colour <- "grey45"
  map_colours  <- grDevices::hcl.colors(120, palette = "Viridis")
  
  sample_50_0 <- NULL
  sample_300_200 <- NULL
  if (show_sample_points) {
    reps_50 <- parameter_df$replicate[parameter_df$config_id == config_50_0 &
                                        parameter_df$converged]
    reps_300 <- parameter_df$replicate[parameter_df$config_id == config_300_200 &
                                         parameter_df$converged]
    if (length(reps_50))
      sample_50_0 <- subsamples[[config_50_0]][[reps_50[1]]]
    if (length(reps_300))
      sample_300_200 <- subsamples[[config_300_200]][[reps_300[1]]]
  }
  
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
    graphics::rect(xmin, ymin, xmax, ymax, col = sea_colour, border = NA)
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
    graphics::text(x,
                   y,
                   label,
                   adj = c(0, 1),
                   font = 2,
                   cex = 1.4)
  }
  
  add_north_arrow <- function() {
    # Bottom-right of panel a, inset slightly to the left.
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
    graphics::text(x,
                   y1 + 0.032 * map_height,
                   "N",
                   font = 2,
                   cex = 0.95)
    graphics::arrows(x, y0, x, y1, length = 0.10, lwd = 1.6)
  }
  
  add_scale_bar <- function() {
    # Fixed 5-km bar. This landscape is projected in metres.
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
  
  draw_panel <- function(r,
                         alpha_r,
                         label,
                         pts = NULL,
                         point_cex = 0.3,
                         cartography = FALSE) {
    graphics::par(
      mar = c(0.05, 0.05, 0.05, 0.05),
      xaxs = "i",
      yaxs = "i"
    )
    draw_basemap()
    
    terra::plot(
      r,
      add = TRUE,
      range = c(0, 1),
      fill_range = TRUE,
      col = map_colours,
      alpha = alpha_r,
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
    
    if (!is.null(pts)) {
      graphics::points(
        pts$x,
        pts$y,
        pch = 21,
        bg = grDevices::adjustcolor("white", alpha.f = 0.90),
        col = "black",
        cex = point_cex,
        lwd = 0.45
      )
    }
    
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
    
    left <- 0.35
    right <- 0.65
    bottom <- 0.48
    top <- 0.70
    breaks <- seq(0, 1, length.out = length(map_colours) + 1)
    
    for (i in seq_along(map_colours)) {
      x1 <- left + (right - left) * breaks[i]
      x2 <- left + (right - left) * breaks[i + 1]
      graphics::rect(x1, bottom, x2, top, col = map_colours[i], border = NA)
    }
    
    graphics::rect(left,
                   bottom,
                   right,
                   top,
                   border = "grey20",
                   lwd = 0.8)
    ticks <- seq(0, 1, by = 0.2)
    tx <- left + (right - left) * ticks
    graphics::segments(tx, bottom, tx, bottom - 0.06, lwd = 0.7)
    graphics::text(
      tx,
      bottom - 0.10,
      labels = sprintf("%.1f", ticks),
      cex = 0.72,
      adj = c(0.5, 1)
    )
    graphics::text(0.5, 0.91, "Occupancy probability", cex = 0.80)
  }
  
  draw_occupancy_figure <- function() {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    
    graphics::layout(
      matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE),
      widths = c(1, 1, 1),
      heights = c(1, 0.11)
    )
    
    draw_panel(true_draw, true_alpha, "a", cartography = TRUE)
    draw_panel(pred_50_draw, pred_50_alpha, "b", sample_50_0, point_cex = 0.43)
    draw_panel(pred_300_draw,
               pred_300_alpha,
               "c",
               sample_300_200,
               point_cex = 0.23)
    draw_colour_bar()
  }
  
  grDevices::png(
    filename = file.path(output_dir, figure_file),
    width = 4800,
    height = 1500,
    units = "px",
    res = 400,
    bg = "white"
  )
  draw_occupancy_figure()
  grDevices::dev.off()
  message("Saved figure: ", file.path(output_dir, figure_file))
  
  if (preview_figure && interactive()) {
    grDevices::dev.new(width = 13, height = 4.1)
    draw_occupancy_figure()
  }
}


# ==============================================================================
# 15. REPRODUCIBILITY INFORMATION
# ==============================================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

# ==============================================================================
# END
# ==============================================================================
