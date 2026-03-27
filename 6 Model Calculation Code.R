# ============================================================
# Script: Model Evaluation and Comparison for 6 Spatial Models
#
# Description:
#   This script trains and evaluates 6 models (LM, RF, BCS, GOS, PS, PSL) 
#   using spatial cross-validation.
#
# Features (Based on manuscript architecture):
#   - LM, RF, BCS, GOS: Uses 10 original covariates.
#   - PS: Uses 40D expanded features with similarity-weighted mean.
#   - PSL (Two-Stream Architecture): 
#       1. Similarity / weights computed on 10D original covariates (h=1.0).
#       2. Local RF trained on 40D expanded features (ntrees=100, mtry=6).
#   - Reproducibility: Global seed set to 42 for RF and PSL.
#   - Performance: Parallel computing enabled for speedup.
#
# Input:
#   data.csv (Requires coordinates, foldID, response variable, and covariates)
#
# Outputs:
#   6models_result.csv  (Predictions and residuals for all models)
#   6models_metrics.csv (R2, RMSE, MAE for each fold and overall)
# ============================================================

# -------------------------
# 0) Load Packages
# -------------------------
# Uncomment the line below to install required packages if not already installed
# install.packages(c("readr", "dplyr", "tibble", "geosimilarity", "ranger"))

library(readr)
library(dplyr)
library(tibble)
library(geosimilarity)
library(ranger)
library(parallel)

# Detect available CPU cores and reserve 1 for the OS
n_cores_available <- max(1L, parallel::detectCores() - 1L)
cat("==================================================\n")
cat("System detected CPU cores. Utilizing", n_cores_available, "cores for parallel computing!\n")
cat("==================================================\n")

# Set global seed for reproducibility
set.seed(42)

# -------------------------
# 1) File Paths
# -------------------------
infile_data    <- "data.csv"
outfile_result <- "6models_result.csv"
outfile_metric <- "6models_metrics.csv"

# -------------------------
# 2) Variable Settings
# -------------------------
covars_10 <- c(
  "Slope", "Biomass", "Precipitation", "TPI", "Temperature",
  "NPP", "SolarRadiation", "WindSpeed", "SinAspect", "CosAspect"
)

features_40 <- c(
  covars_10,
  paste0("C_",  covars_10),
  paste0("IP_", covars_10),
  paste0("IN_", covars_10)
)

response <- "SOC"

# -------------------------
# 3) Model Parameters
# -------------------------
# BCS / GOS
kappa_bcs <- 1.0
kappa_gos <- 0.10

# RF (Global)
rf_trees <- 500L

# PS (40D weighted mean)
ps_params <- list(
  top_k      = 12L,
  h_pattern  = 3,     
  batch_size = 500L
)

# PSL (Two-Stream Architecture: 10D for Similarity, 40D for RF)
psl_params <- list(
  response   = "SOC",
  sim_cols   = covars_10,      # Stream A: 10D original variables for similarity
  rf_cols    = features_40,    # Stream B: 40D expanded features for Local RF
  pool_pct   = 0.30,           # Retain 30% candidate pool
  n_trees    = 100L,           # Number of trees
  h_x        = 1.0,            # Bandwidth parameter (h = 1.0 after standardization)
  mtry_val   = 6L,             # mtry = 6 (floor(sqrt(40)))
  n_cores    = n_cores_available, 
  batch_size = 500L
)

# -------------------------
# 4) Helper: Evaluation Metrics
# -------------------------
calc_metrics_cm <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 2L) {
    return(list(R2 = NA_real_, RMSE = NA_real_, MAE = NA_real_))
  }
  o <- obs[ok]
  p <- pred[ok]
  list(
    R2   = round(1 - sum((o - p)^2) / sum((o - mean(o))^2), 4),
    RMSE = round(sqrt(mean((o - p)^2)), 4),
    MAE  = round(mean(abs(o - p)), 4)
  )
}

# -------------------------
# 5) Helper: Product-Gaussian Similarity
# -------------------------
calc_gos_sim_matrix <- function(train_x, test_x, h = 1, batch_size = 500L) {
  n_tr <- nrow(train_x)
  n_te <- nrow(test_x)
  result <- matrix(0.0, nrow = n_tr, ncol = n_te)
  
  for (b_start in seq(1, n_te, by = batch_size)) {
    b_end  <- min(b_start + batch_size - 1L, n_te)
    b_test <- test_x[b_start:b_end, , drop = FALSE]
    b_size <- b_end - b_start + 1L
    
    d2 <- matrix(0.0, nrow = n_tr, ncol = b_size)
    for (k in seq_len(ncol(train_x))) {
      d2 <- d2 + outer(
        train_x[, k], b_test[, k],
        FUN = function(a, b) ((a - b) / h)^2
      )
    }
    result[, b_start:b_end] <- exp(-0.5 * d2)
  }
  
  result
}

# -------------------------
# 6) PS Model: 40D Similarity-Weighted Mean
# -------------------------
ps_cv <- function(dat, feature_cols, params, response = "SOC") {
  top_k      <- if (!is.null(params$top_k))      params$top_k      else 12L
  h_pattern  <- if (!is.null(params$h_pattern))  params$h_pattern  else 1
  batch_size <- if (!is.null(params$batch_size)) params$batch_size else 300L
  
  folds <- sort(unique(dat$foldID))
  all_results <- vector("list", length(folds))
  
  for (fi in seq_along(folds)) {
    fold <- folds[fi]
    cat(sprintf("PS: fold %d / %d ...\n", fi, length(folds)))
    
    train_idx <- which(dat$foldID != fold)
    test_idx  <- which(dat$foldID == fold)
    
    train_df <- dat[train_idx, , drop = FALSE]
    test_df  <- dat[test_idx,  , drop = FALSE]
    
    tr_means <- colMeans(train_df[, feature_cols, drop = FALSE], na.rm = TRUE)
    tr_sds   <- apply(train_df[, feature_cols, drop = FALSE], 2L, sd, na.rm = TRUE)
    tr_sds[!is.finite(tr_sds) | tr_sds == 0] <- 1
    
    train_sc <- scale(train_df[, feature_cols, drop = FALSE], center = tr_means, scale = tr_sds)
    test_sc  <- scale(test_df[, feature_cols, drop = FALSE],  center = tr_means, scale = tr_sds)
    
    sim_mat <- calc_gos_sim_matrix(
      train_x    = train_sc,
      test_x     = test_sc,
      h          = h_pattern,
      batch_size = batch_size
    )
    
    n_train <- nrow(train_df)
    k_use   <- min(top_k, n_train)
    
    preds <- vapply(seq_len(ncol(sim_mat)), function(j) {
      s   <- sim_mat[, j]
      ord <- order(s, decreasing = TRUE)[seq_len(k_use)]
      w   <- s[ord]
      y   <- train_df[[response]][ord]
      
      if (!all(is.finite(w)) || sum(w, na.rm = TRUE) <= 0) {
        return(mean(y, na.rm = TRUE))
      }
      sum(w * y, na.rm = TRUE) / sum(w, na.rm = TRUE)
    }, numeric(1))
    
    all_results[[fi]] <- data.frame(
      foldID   = fold,
      point_id = test_idx,
      obs      = test_df[[response]],
      pred_ps  = preds,
      stringsAsFactors = FALSE
    )
  }
  
  result <- do.call(rbind, all_results)
  rownames(result) <- NULL
  
  rows <- list()
  for (f in folds) {
    sub <- result[result$foldID == f, ]
    m <- calc_metrics_cm(sub$obs, sub$pred_ps)
    rows[[length(rows) + 1L]] <- data.frame(
      Fold = as.character(f), Model = "PS",
      R2 = m$R2, RMSE = m$RMSE, MAE = m$MAE,
      stringsAsFactors = FALSE
    )
  }
  m_all <- calc_metrics_cm(result$obs, result$pred_ps)
  rows[[length(rows) + 1L]] <- data.frame(
    Fold = "Overall", Model = "PS",
    R2 = m_all$R2, RMSE = m_all$RMSE, MAE = m_all$MAE,
    stringsAsFactors = FALSE
  )
  
  attr(result, "fold_metrics") <- do.call(rbind, rows)
  result
}

# -------------------------
# 7) PSL Model: Two-Stream Architecture
# -------------------------
psl_cv <- function(dat, params) {
  response   <- params$response
  sim_cols   <- params$sim_cols   
  rf_cols    <- params$rf_cols    
  pool_pct   <- params$pool_pct   
  n_trees    <- params$n_trees    
  h_x        <- params$h_x        
  mtry_val   <- params$mtry_val   
  n_cores    <- params$n_cores
  batch_size <- params$batch_size
  
  folds <- sort(unique(dat$foldID))
  all_results <- vector("list", length(folds))
  
  for (fi in seq_along(folds)) {
    fold <- folds[fi]
    cat(sprintf("PSL (Two-Stream): fold %d / %d ...\n", fi, length(folds)))
    
    train_idx <- which(dat$foldID != fold)
    test_idx  <- which(dat$foldID == fold)
    
    train_df <- dat[train_idx, , drop = FALSE]
    test_df  <- dat[test_idx,  , drop = FALSE]
    
    n_train <- nrow(train_df)
    n_test  <- nrow(test_df)
    
    # Core 1: Standardize in 10D space and compute geographical similarity
    tr_means <- colMeans(train_df[, sim_cols, drop = FALSE], na.rm = TRUE)
    tr_sds   <- apply(train_df[, sim_cols, drop = FALSE], 2L, sd, na.rm = TRUE)
    tr_sds[!is.finite(tr_sds) | tr_sds == 0] <- 1
    
    train_sc <- scale(train_df[, sim_cols, drop = FALSE], center = tr_means, scale = tr_sds)
    test_sc  <- scale(test_df[, sim_cols, drop = FALSE],  center = tr_means, scale = tr_sds)
    
    sim_mat <- calc_gos_sim_matrix(
      train_x    = train_sc,
      test_x     = test_sc,
      h          = h_x,
      batch_size = batch_size
    )
    
    n_pool <- max(10L, floor(n_train * pool_pct))
    
    # Core 2: Extract features in 40D space for local training
    train_feat <- train_df[, c(rf_cols, response), drop = FALSE]
    test_feat  <- test_df[, rf_cols, drop = FALSE]
    
    rf_formula <- as.formula(paste(response, "~ ."))
    
    preds <- mclapply(seq_len(n_test), function(j) {
      sim_j <- sim_mat[, j]
      
      # Select Top pool based on 10D similarity
      pool_idx <- order(sim_j, decreasing = TRUE)[seq_len(n_pool)]
      local_df <- train_feat[pool_idx, , drop = FALSE]
      
      if (nrow(local_df) < 5L) return(NA_real_)
      
      # Utilize 10D similarity weights directly
      w <- sim_j[pool_idx] 
      
      fit <- ranger(
        formula       = rf_formula,
        data          = local_df,
        case.weights  = w,          
        num.trees     = n_trees,    
        mtry          = mtry_val,   
        min.node.size = 5L,
        num.threads   = 1L,         # Must be 1 inside parallel loop
        verbose       = FALSE,
        seed          = 42          # Aligned seed
      )
      
      # Predict using the 40D features of the target point
      predict(fit, data = test_feat[j, , drop = FALSE])$predictions
    }, mc.cores = n_cores)      
    
    all_results[[fi]] <- data.frame(
      foldID   = fold,
      point_id = test_idx,
      obs      = test_df[[response]],
      pred_psl = unlist(preds),
      stringsAsFactors = FALSE
    )
  }
  
  result <- do.call(rbind, all_results)
  rownames(result) <- NULL
  
  rows <- list()
  for (f in folds) {
    sub <- result[result$foldID == f, ]
    m <- calc_metrics_cm(sub$obs, sub$pred_psl)
    rows[[length(rows) + 1L]] <- data.frame(
      Fold = as.character(f), Model = "PSL",
      R2 = m$R2, RMSE = m$RMSE, MAE = m$MAE,
      stringsAsFactors = FALSE
    )
  }
  m_all <- calc_metrics_cm(result$obs, result$pred_psl)
  rows[[length(rows) + 1L]] <- data.frame(
    Fold = "Overall", Model = "PSL",
    R2 = m_all$R2, RMSE = m_all$RMSE, MAE = m_all$MAE,
    stringsAsFactors = FALSE
  )
  
  attr(result, "fold_metrics") <- do.call(rbind, rows)
  result
}

# -------------------------
# 8) Read Data
# -------------------------
need_cols <- c(
  "longitude", "latitude", "hex_id", "foldID", response,
  covars_10, features_40
)

dat <- read_csv(infile_data, show_col_types = FALSE)

miss_cols <- setdiff(unique(need_cols), names(dat))
if (length(miss_cols) > 0) {
  stop("Missing columns in data file: ", paste(miss_cols, collapse = ", "))
}

dat <- dat %>%
  filter(if_all(all_of(c(response, covars_10, features_40)), is.finite)) %>%
  filter(!is.na(foldID), !is.na(hex_id)) %>%
  mutate(
    point_id = row_number(),
    foldID   = as.integer(foldID)
  )

folds <- sort(unique(dat$foldID))

cat("Loaded ", nrow(dat), " rows across ", length(folds), " folds.\n", sep = "")
cat("Configuration Summary:\n")
cat(" - LM/RF/BCS/GOS use 10 original covariates.\n")
cat(" - PS uses 40 expanded features.\n")
cat(" - PSL uses 10D for similarity weighting, and 40D for Local RF training.\n\n")

# -------------------------
# 9) Run LM / RF / BCS / GOS on 10 Covariates
# -------------------------
results_10d <- vector("list", length(folds))

for (fi in seq_along(folds)) {
  f <- folds[fi]
  cat("Running 10D models: fold ", f, " / ", length(folds), " ...\n", sep = "")
  
  train_idx <- which(dat$foldID != f)
  test_idx  <- which(dat$foldID == f)
  
  train_df <- dat[train_idx, c(response, covars_10)] %>% as.data.frame()
  test_df  <- dat[test_idx,  covars_10] %>% as.data.frame()
  
  # BCS
  bcs_res <- gos(
    SOC ~ .,
    data    = train_df,
    newdata = test_df,
    kappa   = kappa_bcs,
    cores   = n_cores_available
  )
  pred_bcs <- as.numeric(bcs_res$pred)
  
  # GOS
  gos_res <- gos(
    SOC ~ .,
    data    = train_df,
    newdata = test_df,
    kappa   = kappa_gos,
    cores   = n_cores_available
  )
  pred_gos <- as.numeric(gos_res$pred)
  
  # LM
  lm_fit  <- lm(SOC ~ ., data = train_df)
  pred_lm <- as.numeric(predict(lm_fit, newdata = test_df))
  
  # RF (Global)
  rf_fit <- ranger(
    SOC ~ .,
    data        = train_df,
    num.trees   = rf_trees,
    seed        = 42,
    num.threads = n_cores_available
  )
  pred_rf <- as.numeric(predict(rf_fit, data = test_df)$predictions)
  
  results_10d[[fi]] <- tibble(
    point_id  = dat$point_id[test_idx],
    foldID    = f,
    obs       = dat[[response]][test_idx],
    pred_bcs  = pred_bcs,
    pred_gos  = pred_gos,
    pred_lm   = pred_lm,
    pred_rf   = pred_rf
  )
}

cv_10d <- bind_rows(results_10d)

# -------------------------
# 10) Run PS & PSL
# -------------------------
cat("\nRunning PS on 40D features ...\n")
ps_res <- ps_cv(
  dat         = dat,
  feature_cols = features_40,
  params      = ps_params,
  response    = response
)

cat("\nRunning PSL with Two-Stream Architecture ...\n")
psl_res <- psl_cv(
  dat         = dat,
  params      = psl_params
)

# -------------------------
# 11) Merge All Predictions
# -------------------------
res_all <- dat %>%
  select(point_id, foldID, hex_id, longitude, latitude, obs = all_of(response)) %>%
  left_join(cv_10d,  by = c("point_id", "foldID", "obs")) %>%
  left_join(ps_res  %>% select(point_id, pred_ps),  by = "point_id") %>%
  left_join(psl_res %>% select(point_id, pred_psl), by = "point_id") %>%
  mutate(
    resid_bcs = obs - pred_bcs,
    resid_gos = obs - pred_gos,
    resid_lm  = obs - pred_lm,
    resid_rf  = obs - pred_rf,
    resid_ps  = obs - pred_ps,
    resid_psl = obs - pred_psl
  ) %>%
  arrange(point_id)

# -------------------------
# 12) Save Results
# -------------------------
write_csv(res_all, outfile_result)
cat("\nSaved result file:\n", outfile_result, "\n", sep = "")

# -------------------------
# 13) Build Metrics Summary
# -------------------------
model_cols <- c(
  BCS = "pred_bcs",
  GOS = "pred_gos",
  LM  = "pred_lm",
  RF  = "pred_rf",
  PS  = "pred_ps",
  PSL = "pred_psl"
)

rows <- list()

for (f in folds) {
  sub <- res_all %>% filter(foldID == f)
  for (mn in names(model_cols)) {
    m <- calc_metrics_cm(sub$obs, sub[[model_cols[mn]]])
    rows[[length(rows) + 1L]] <- tibble(
      Fold  = as.character(f),
      Model = mn,
      R2    = m$R2,
      RMSE  = m$RMSE,
      MAE   = m$MAE
    )
  }
}

for (mn in names(model_cols)) {
  m <- calc_metrics_cm(res_all$obs, res_all[[model_cols[mn]]])
  rows[[length(rows) + 1L]] <- tibble(
    Fold  = "Overall",
    Model = mn,
    R2    = m$R2,
    RMSE  = m$RMSE,
    MAE   = m$MAE
  )
}

metrics_all <- bind_rows(rows)

write_csv(metrics_all, outfile_metric)
cat("Saved metrics file:\n", outfile_metric, "\n", sep = "")

cat("\n--- Overall Pooled Metrics ---\n")
print(metrics_all %>% filter(Fold == "Overall"))