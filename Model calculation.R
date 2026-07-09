# ============================================================
# Build 40D features and evaluate six spatial models:
# SGWR, GOS, LM, RF, PS, and PSL
#
# This script corresponds to the final PSL algorithm used in the manuscript.
#
# PSL logic:
#   Stream A:
#     Use the original 10-dimensional covariate space to compute
#     product-Gaussian similarity and select the most similar training samples.
#
#   Stream B:
#     Use the expanded 40-dimensional pattern-enhanced feature space
#     to train a similarity-weighted local Random Forest.
#
# Required input:
#   data.csv
#
# Required columns:
#   longitude, latitude, SOC, hex_id, foldID,
#   and the 10 environmental covariates listed below.
#
# Output:
#   6models_residuals.csv
# ============================================================

# -------------------------
# 0) Load packages
# -------------------------
library(readr)
library(dplyr)
library(tibble)
library(sf)
library(spdep)
library(geocomplexity)
library(geosimilarity)
library(ranger)
library(parallel)

# -------------------------
# 1) Input and output files
# -------------------------
# Place data.csv in the same folder as this script.
infile_data <- "data.csv"
outfile_residual <- "6models_residuals.csv"

# -------------------------
# 2) Variables
# -------------------------
response <- "SOC"

# Original 10-dimensional environmental covariates.
covars_10 <- c(
  "Slope", "Biomass", "Precipitation", "TPI", "Temperature",
  "NPP", "SolarRadiation", "WindSpeed", "SinAspect", "CosAspect"
)

# Expanded 40-dimensional feature space:
# 10 original covariates + 10 geocomplexity features +
# 10 positive outlier-strength features + 10 negative outlier-strength features.
features_40 <- c(
  covars_10,
  paste0("C_",  covars_10),
  paste0("IP_", covars_10),
  paste0("IN_", covars_10)
)

need_cols_raw <- c(
  "longitude", "latitude", response,
  covars_10,
  "hex_id", "foldID"
)

# -------------------------
# 3) General parameters
# -------------------------
set.seed(42)

n_cores_available <- max(1L, parallel::detectCores() - 1L)
cat("Using", n_cores_available, "CPU cores.\n")

# Neighborhood size used to construct spatial pattern features.
k_nn <- 12L

# GOS similarity parameter.
kappa_gos <- 0.10

# SGWR settings:
# geo_k: local geographic neighborhood size
# h_attr: attribute similarity bandwidth
# ridge: ridge stabilization for local regression
# min_n: minimum local sample size
sgwr_params <- list(
  geo_k = 300L,
  h_attr = 1.0,
  ridge = 1e-6,
  min_n = 30L
)

# Tuned global RF baseline.
rf_params <- list(
  num_trees = 1500L,
  mtry = 2L,
  min_node_size = 1L
)

# PS uses similarity-weighted averaging in the 40D feature space.
ps_params <- list(
  top_k = 12L,
  h_pattern = 3,
  batch_size = 500L
)

# PSL uses 10D covariates for similarity search and 40D features for local RF learning.
# kappa_psl is the proportion of most similar training samples retained for local RF.
psl_params <- list(
  response = response,
  sim_cols = covars_10,
  rf_cols = features_40,
  kappa_psl = 0.70,
  n_trees = 100L,
  h_x = 1.0,
  mtry_val = 6L,
  min_node_size = 5L,
  n_cores = n_cores_available,
  batch_size = 500L
)

# -------------------------
# 4) Read and clean input data
# -------------------------
dat_raw <- read_csv(infile_data, show_col_types = FALSE)

miss_cols <- setdiff(need_cols_raw, names(dat_raw))
if (length(miss_cols) > 0) {
  stop("Missing columns in data file: ", paste(miss_cols, collapse = ", "))
}

dat_raw <- dat_raw %>%
  filter(if_all(all_of(c("longitude", "latitude", response, covars_10)), is.finite)) %>%
  filter(!is.na(hex_id), !is.na(foldID)) %>%
  mutate(
    point_id = row_number(),
    foldID = as.integer(foldID)
  )

n <- nrow(dat_raw)

if (n < k_nn + 1L) {
  stop("Too few rows for k_nn = ", k_nn, ".")
}

cat("Loaded ", n, " rows.\n", sep = "")

# Coordinates are assumed to be projected coordinates in EPSG:5070.
# If longitude and latitude are geographic coordinates, they should be projected
# before running this script.
pts <- st_as_sf(
  dat_raw,
  coords = c("longitude", "latitude"),
  crs = 5070,
  remove = FALSE
)

xy <- st_coordinates(pts)

# -------------------------
# 5) Build 40D spatial pattern features
# -------------------------

# ---------- 5.1 Geocomplexity features ----------
nb_c <- knn2nb(knearneigh(xy, k = k_nn))
wt_c <- nb2mat(nb_c, style = "B", zero.policy = TRUE)

calc_gc_one <- function(varname, sf_points, wt_mat) {
  sf_tmp <- sf_points[, varname, drop = FALSE]
  
  gc_res <- geocd_vector(
    sfj = sf_tmp,
    wt = wt_mat,
    method = "moran",
    normalize = TRUE,
    returnsf = FALSE
  )
  
  gc_res <- as.data.frame(gc_res)
  gc_col <- grep("^GC", names(gc_res), value = TRUE)
  
  if (length(gc_col) == 0) {
    stop("No GC column returned for variable: ", varname)
  }
  
  as.numeric(gc_res[[gc_col[1]]])
}

cat("Computing geocomplexity features ...\n")

C_mat <- sapply(covars_10, function(v) calc_gc_one(v, pts, wt_c))
C_mat <- as.data.frame(C_mat)
names(C_mat) <- paste0("C_", covars_10)

cat("Computed C_* features.\n")

# ---------- 5.2 Positive and negative local outlier-strength features ----------
cat("Computing positive and negative outlier-strength features ...\n")

knn_idx <- knearneigh(xy, k = k_nn)$nn

IP_mat <- matrix(0.0, nrow = n, ncol = length(covars_10))
IN_mat <- matrix(0.0, nrow = n, ncol = length(covars_10))

colnames(IP_mat) <- paste0("IP_", covars_10)
colnames(IN_mat) <- paste0("IN_", covars_10)

for (ki in seq_along(covars_10)) {
  vname <- covars_10[ki]
  xk <- dat_raw[[vname]]
  
  nn_mat <- matrix(xk[knn_idx], nrow = n, ncol = k_nn)
  
  mu_k <- rowMeans(nn_mat, na.rm = TRUE)
  dev <- sweep(nn_mat, 1L, mu_k, "-")
  sd_k <- sqrt(rowSums(dev^2, na.rm = TRUE) / (k_nn - 1L))
  
  IP_mat[, ki] <- pmax(0, xk - mu_k - 2 * sd_k)
  IN_mat[, ki] <- pmax(0, mu_k - 2 * sd_k - xk)
}

IP_mat <- as.data.frame(IP_mat)
IN_mat <- as.data.frame(IN_mat)

cat("Computed IP_* and IN_* features.\n")

# ---------- 5.3 Combine original variables and spatial pattern features ----------
dat_40d <- dat_raw %>%
  select(
    point_id,
    longitude, latitude, all_of(response),
    all_of(covars_10),
    hex_id, foldID
  ) %>%
  bind_cols(C_mat, IP_mat, IN_mat)

if (!all(features_40 %in% names(dat_40d))) {
  stop("Some 40D feature columns are missing.")
}

cat("40D feature block verified.\n")

# -------------------------
# 6) Evaluation metrics
# -------------------------
calc_metrics_cm <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  
  if (sum(ok) < 2L) {
    return(list(R2 = NA_real_, RMSE = NA_real_, MAE = NA_real_))
  }
  
  o <- obs[ok]
  p <- pred[ok]
  
  list(
    R2 = round(1 - sum((o - p)^2) / sum((o - mean(o))^2), 4),
    RMSE = round(sqrt(mean((o - p)^2)), 4),
    MAE = round(mean(abs(o - p)), 4)
  )
}

# -------------------------
# 7) Product-Gaussian similarity matrix
# -------------------------
calc_gos_sim_matrix <- function(train_x, test_x, h = 1.0, batch_size = 500L) {
  n_tr <- nrow(train_x)
  n_te <- nrow(test_x)
  
  result <- matrix(0.0, nrow = n_tr, ncol = n_te)
  
  for (b_start in seq(1, n_te, by = batch_size)) {
    b_end <- min(b_start + batch_size - 1L, n_te)
    b_test <- test_x[b_start:b_end, , drop = FALSE]
    b_size <- b_end - b_start + 1L
    
    d2 <- matrix(0.0, nrow = n_tr, ncol = b_size)
    
    for (k in seq_len(ncol(train_x))) {
      d2 <- d2 + outer(
        train_x[, k],
        b_test[, k],
        FUN = function(a, b) ((a - b) / h)^2
      )
    }
    
    result[, b_start:b_end] <- exp(-0.5 * d2)
  }
  
  result
}

# -------------------------
# 8) SGWR prediction for one fold
# -------------------------
sgwr_predict_fold <- function(train_df, test_df, covars, response, params) {
  geo_k <- params$geo_k
  h_attr <- params$h_attr
  ridge <- params$ridge
  min_n <- params$min_n
  
  n_train <- nrow(train_df)
  n_test <- nrow(test_df)
  
  geo_k <- min(geo_k, n_train)
  
  tr_means <- colMeans(train_df[, covars, drop = FALSE], na.rm = TRUE)
  tr_sds <- apply(train_df[, covars, drop = FALSE], 2L, sd, na.rm = TRUE)
  tr_sds[!is.finite(tr_sds) | tr_sds == 0] <- 1
  
  train_x <- scale(train_df[, covars, drop = FALSE], center = tr_means, scale = tr_sds)
  test_x <- scale(test_df[, covars, drop = FALSE], center = tr_means, scale = tr_sds)
  
  train_xy <- as.matrix(train_df[, c("longitude", "latitude")])
  test_xy <- as.matrix(test_df[, c("longitude", "latitude")])
  
  y_train <- train_df[[response]]
  
  preds <- numeric(n_test)
  
  for (j in seq_len(n_test)) {
    d_geo <- sqrt(
      (train_xy[, 1] - test_xy[j, 1])^2 +
        (train_xy[, 2] - test_xy[j, 2])^2
    )
    
    geo_idx <- order(d_geo)[seq_len(geo_k)]
    
    d_local <- d_geo[geo_idx]
    bw_geo <- max(d_local, na.rm = TRUE)
    
    if (!is.finite(bw_geo) || bw_geo == 0) {
      bw_geo <- 1
    }
    
    w_geo <- exp(-0.5 * (d_local / bw_geo)^2)
    
    dx_attr <- sweep(
      train_x[geo_idx, , drop = FALSE],
      2L,
      test_x[j, ],
      "-"
    )
    
    d2_attr <- rowSums(dx_attr^2)
    w_attr <- exp(-0.5 * d2_attr / (h_attr^2))
    
    w <- w_geo * w_attr
    
    ok <- is.finite(w) & is.finite(y_train[geo_idx]) & w > 0
    
    if (sum(ok) < min_n) {
      preds[j] <- mean(y_train[geo_idx], na.rm = TRUE)
      next
    }
    
    use_idx <- geo_idx[ok]
    w_use <- w[ok]
    
    local_df <- train_df[use_idx, c(response, covars), drop = FALSE]
    
    X <- model.matrix(
      reformulate(covars),
      data = local_df
    )
    y <- local_df[[response]]
    
    W <- sqrt(w_use)
    Xw <- X * W
    yw <- y * W
    
    XtX <- crossprod(Xw)
    Xty <- crossprod(Xw, yw)
    
    ridge_mat <- diag(ridge, ncol(XtX))
    ridge_mat[1, 1] <- 0
    
    beta <- tryCatch(
      solve(XtX + ridge_mat, Xty),
      error = function(e) NULL
    )
    
    if (is.null(beta) || any(!is.finite(beta))) {
      preds[j] <- weighted.mean(y, w_use, na.rm = TRUE)
    } else {
      x0 <- model.matrix(
        reformulate(covars),
        data = test_df[j, covars, drop = FALSE]
      )
      preds[j] <- as.numeric(x0 %*% beta)
    }
  }
  
  preds
}

# -------------------------
# 9) PS model
# -------------------------
# PS uses similarity-weighted averaging based on the expanded 40D feature space.
ps_cv <- function(dat, feature_cols, params, response = "SOC") {
  top_k <- params$top_k
  h_pattern <- params$h_pattern
  batch_size <- params$batch_size
  
  folds <- sort(unique(dat$foldID))
  all_results <- vector("list", length(folds))
  
  for (fi in seq_along(folds)) {
    fold <- folds[fi]
    cat(sprintf("PS: fold %d / %d ...\n", fi, length(folds)))
    
    train_idx <- which(dat$foldID != fold)
    test_idx <- which(dat$foldID == fold)
    
    train_df <- dat[train_idx, , drop = FALSE]
    test_df <- dat[test_idx, , drop = FALSE]
    
    tr_means <- colMeans(train_df[, feature_cols, drop = FALSE], na.rm = TRUE)
    tr_sds <- apply(train_df[, feature_cols, drop = FALSE], 2L, sd, na.rm = TRUE)
    tr_sds[!is.finite(tr_sds) | tr_sds == 0] <- 1
    
    train_sc <- scale(
      train_df[, feature_cols, drop = FALSE],
      center = tr_means,
      scale = tr_sds
    )
    
    test_sc <- scale(
      test_df[, feature_cols, drop = FALSE],
      center = tr_means,
      scale = tr_sds
    )
    
    sim_mat <- calc_gos_sim_matrix(
      train_x = train_sc,
      test_x = test_sc,
      h = h_pattern,
      batch_size = batch_size
    )
    
    n_train <- nrow(train_df)
    k_use <- min(top_k, n_train)
    
    preds <- vapply(seq_len(ncol(sim_mat)), function(j) {
      s <- sim_mat[, j]
      ord <- order(s, decreasing = TRUE)[seq_len(k_use)]
      w <- s[ord]
      y <- train_df[[response]][ord]
      
      if (!all(is.finite(w)) || sum(w, na.rm = TRUE) <= 0) {
        return(mean(y, na.rm = TRUE))
      }
      
      sum(w * y, na.rm = TRUE) / sum(w, na.rm = TRUE)
    }, numeric(1))
    
    all_results[[fi]] <- data.frame(
      foldID = fold,
      point_id = test_df$point_id,
      obs = test_df[[response]],
      pred_ps = preds,
      stringsAsFactors = FALSE
    )
  }
  
  result <- do.call(rbind, all_results)
  rownames(result) <- NULL
  result
}

# -------------------------
# 10) PSL model
# -------------------------
# PSL uses 10D covariates for similarity search and 40D features for local RF learning.
psl_cv <- function(dat, params) {
  response <- params$response
  sim_cols <- params$sim_cols
  rf_cols <- params$rf_cols
  kappa_psl <- params$kappa_psl
  n_trees <- params$n_trees
  h_x <- params$h_x
  mtry_val <- params$mtry_val
  min_node_size <- params$min_node_size
  n_cores <- params$n_cores
  batch_size <- params$batch_size
  
  folds <- sort(unique(dat$foldID))
  all_results <- vector("list", length(folds))
  
  for (fi in seq_along(folds)) {
    fold <- folds[fi]
    cat(sprintf("PSL: fold %d / %d ...\n", fi, length(folds)))
    
    train_idx <- which(dat$foldID != fold)
    test_idx <- which(dat$foldID == fold)
    
    train_df <- dat[train_idx, , drop = FALSE]
    test_df <- dat[test_idx, , drop = FALSE]
    
    n_train <- nrow(train_df)
    n_test <- nrow(test_df)
    
    # Stream A: 10D similarity search.
    tr_means <- colMeans(train_df[, sim_cols, drop = FALSE], na.rm = TRUE)
    tr_sds <- apply(train_df[, sim_cols, drop = FALSE], 2L, sd, na.rm = TRUE)
    tr_sds[!is.finite(tr_sds) | tr_sds == 0] <- 1
    
    train_sc <- scale(
      train_df[, sim_cols, drop = FALSE],
      center = tr_means,
      scale = tr_sds
    )
    
    test_sc <- scale(
      test_df[, sim_cols, drop = FALSE],
      center = tr_means,
      scale = tr_sds
    )
    
    sim_mat <- calc_gos_sim_matrix(
      train_x = train_sc,
      test_x = test_sc,
      h = h_x,
      batch_size = batch_size
    )
    
    n_pool <- max(10L, floor(n_train * kappa_psl))
    n_pool <- min(n_pool, n_train)
    
    train_feat <- train_df[, c(rf_cols, response), drop = FALSE]
    test_feat <- test_df[, rf_cols, drop = FALSE]
    
    rf_formula <- as.formula(paste(response, "~ ."))
    
    # Stream B: similarity-weighted local Random Forest.
    preds <- mclapply(seq_len(n_test), function(j) {
      sim_j <- sim_mat[, j]
      pool_idx <- order(sim_j, decreasing = TRUE)[seq_len(n_pool)]
      
      local_df <- train_feat[pool_idx, , drop = FALSE]
      w <- sim_j[pool_idx]
      
      if (nrow(local_df) < 5L || sum(w, na.rm = TRUE) <= 0) {
        return(mean(local_df[[response]], na.rm = TRUE))
      }
      
      fit <- ranger(
        formula = rf_formula,
        data = local_df,
        case.weights = w,
        num.trees = n_trees,
        mtry = mtry_val,
        min.node.size = min_node_size,
        num.threads = 1L,
        verbose = FALSE,
        seed = 42
      )
      
      as.numeric(
        predict(
          fit,
          data = test_feat[j, , drop = FALSE]
        )$predictions
      )
    }, mc.cores = n_cores)
    
    all_results[[fi]] <- data.frame(
      foldID = fold,
      point_id = test_df$point_id,
      obs = test_df[[response]],
      pred_psl = unlist(preds),
      stringsAsFactors = FALSE
    )
  }
  
  result <- do.call(rbind, all_results)
  rownames(result) <- NULL
  result
}

# -------------------------
# 11) Run LM, RF, SGWR, and GOS
# -------------------------
dat <- dat_40d
folds <- sort(unique(dat$foldID))

cat("Running LM, RF, SGWR, and GOS ...\n")

results_10d <- vector("list", length(folds))

for (fi in seq_along(folds)) {
  f <- folds[fi]
  
  cat("Fold ", fi, " / ", length(folds), " | foldID = ", f, " ...\n", sep = "")
  
  train_idx <- which(dat$foldID != f)
  test_idx <- which(dat$foldID == f)
  
  train_df <- dat[train_idx, c(response, covars_10)] %>%
    as.data.frame()
  
  test_df <- dat[test_idx, covars_10] %>%
    as.data.frame()
  
  sgwr_train_df <- dat[train_idx, c(response, "longitude", "latitude", covars_10)] %>%
    as.data.frame()
  
  sgwr_test_df <- dat[test_idx, c("longitude", "latitude", covars_10)] %>%
    as.data.frame()
  
  pred_sgwr <- sgwr_predict_fold(
    train_df = sgwr_train_df,
    test_df = sgwr_test_df,
    covars = covars_10,
    response = response,
    params = sgwr_params
  )
  
  pred_gos <- as.numeric(
    gos(
      SOC ~ .,
      data = train_df,
      newdata = test_df,
      kappa = kappa_gos,
      cores = n_cores_available
    )$pred
  )
  
  lm_fit <- lm(SOC ~ ., data = train_df)
  pred_lm <- as.numeric(predict(lm_fit, newdata = test_df))
  
  rf_fit <- ranger(
    SOC ~ .,
    data = train_df,
    num.trees = rf_params$num_trees,
    mtry = rf_params$mtry,
    min.node.size = rf_params$min_node_size,
    seed = 42,
    num.threads = n_cores_available
  )
  
  pred_rf <- as.numeric(predict(rf_fit, data = test_df)$predictions)
  
  results_10d[[fi]] <- tibble(
    point_id = dat$point_id[test_idx],
    foldID = f,
    obs = dat[[response]][test_idx],
    pred_sgwr = pred_sgwr,
    pred_gos = pred_gos,
    pred_lm = pred_lm,
    pred_rf = pred_rf
  )
}

cv_10d <- bind_rows(results_10d)

# -------------------------
# 12) Run PS
# -------------------------
cat("Running PS ...\n")

ps_res <- ps_cv(
  dat = dat,
  feature_cols = features_40,
  params = ps_params,
  response = response
)

# -------------------------
# 13) Run PSL
# -------------------------
cat("Running PSL ...\n")

psl_res <- psl_cv(
  dat = dat,
  params = psl_params
)

# -------------------------
# 14) Merge predictions and calculate residuals
# -------------------------
res_all <- dat %>%
  select(point_id, foldID, hex_id, longitude, latitude, obs = all_of(response)) %>%
  left_join(cv_10d, by = c("point_id", "foldID", "obs")) %>%
  left_join(ps_res %>% select(point_id, pred_ps), by = "point_id") %>%
  left_join(psl_res %>% select(point_id, pred_psl), by = "point_id") %>%
  mutate(
    resid_sgwr = obs - pred_sgwr,
    resid_gos = obs - pred_gos,
    resid_lm = obs - pred_lm,
    resid_rf = obs - pred_rf,
    resid_ps = obs - pred_ps,
    resid_psl = obs - pred_psl
  ) %>%
  arrange(point_id)

write_csv(res_all, outfile_residual)

# -------------------------
# 15) Summarize model performance
# -------------------------
model_cols <- c(
  SGWR = "pred_sgwr",
  GOS = "pred_gos",
  LM = "pred_lm",
  RF = "pred_rf",
  PS = "pred_ps",
  PSL = "pred_psl"
)

rows <- list()

for (f in folds) {
  sub <- res_all %>% filter(foldID == f)
  
  for (mn in names(model_cols)) {
    m <- calc_metrics_cm(sub$obs, sub[[model_cols[mn]]])
    
    rows[[length(rows) + 1L]] <- tibble(
      Fold = as.character(f),
      Model = mn,
      R2 = m$R2,
      RMSE = m$RMSE,
      MAE = m$MAE
    )
  }
}

for (mn in names(model_cols)) {
  m <- calc_metrics_cm(res_all$obs, res_all[[model_cols[mn]]])
  
  rows[[length(rows) + 1L]] <- tibble(
    Fold = "Overall",
    Model = mn,
    R2 = m$R2,
    RMSE = m$RMSE,
    MAE = m$MAE
  )
}

metrics_all <- bind_rows(rows)

cat("\n--- Overall Metrics ---\n")
print(metrics_all %>% filter(Fold == "Overall"))

cat("\n--- Fold-level Metrics ---\n")
print(metrics_all %>% filter(Fold != "Overall"))

cat("\nSaved residual file to:\n")
cat(outfile_residual, "\n")

cat("\nDone.\n")