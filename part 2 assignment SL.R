# --------------------------------------------------------------------
# Part 1 — Assignment Statistical Learning
# Authors:
# Group 6:
# Alexandre Lavrinenko    xxxxxx
# Ensar Tasgin:           646820
# Sanne Maasman           xxxxxx 
# Sebastiaan van Helden   822236
# --------------------------------------------------------------------

hybrid_kNN <- function(X, k, ntree = 500, seed = NULL) {
  
  X <- as.matrix(X)
  p <- ncol(X)
  
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(p))
  }
  
  X_original <- X
  X_hat <- X
  missing <- is.na(X)
  weights <- list()
  
  # Robust scaling for kNN distances
  scales <- apply(X, 2, IQR, na.rm = TRUE)
  
  # Use SD if IQR = 0
  zero_iqr <- scales == 0
  scales[zero_iqr] <-
    apply(X[, zero_iqr, drop = FALSE], 2, sd, na.rm = TRUE)
  
  # Constant variables do not need scaling
  scales[is.na(scales) | scales == 0] <- 1
  
  X_scaled <- sweep(X, 2, scales, "/")
  
  
  # Random-forest importance weights
 
  get_rf_weights <- function(target) {
    
    predictors <- setdiff(seq_len(p), target)
    predictor_names <- colnames(X)[predictors]
    
    # Only observations with observed response can train the RF
    rf_rows <- !is.na(X[, target])
    
    y_rf <- X[rf_rows, target]
    x_rf <- X[rf_rows, predictors, drop = FALSE]
    
    
    # Median-impute missing RF predictors
    medians <- apply(
      X[, predictors, drop = FALSE],
      2,
      median,
      na.rm = TRUE
    )
    
    for (j in seq_along(predictors)) {
      x_rf[is.na(x_rf[, j]), j] <- medians[j]
    }
    
    
    if (!is.null(seed)) {
      set.seed(seed + target)
    }
    
    fit <- randomForest::randomForest(
      x = x_rf,
      y = y_rf,
      ntree = ntree,
      importance = TRUE
    )
    
    # Permutation importance
    importance <- drop(
      randomForest::importance(fit, type = 1)
    )
    
    names(importance) <- predictor_names
    
    # Distance weights must be non-negative
    importance <- pmax(importance, 0)
    
    # If RF gives no positive importance, use equal weights
    if (sum(importance) == 0) {
      importance[] <- 1
    }
    
    importance / sum(importance)
  }
  
  
  
  # Hybrid kNN imputation
  
  for (target in seq_len(p)) {
    
    recipients <- which(missing[, target])
    
    if (length(recipients) == 0) {
      next
    }
    
    donors <- which(!missing[, target])
    predictors <- setdiff(seq_len(p), target)
    
    rf_weights <- get_rf_weights(target)
    weights_name <- colnames(X)[target]
    
    
    weights[[weights_name]] <- rf_weights
    
    
    for (recipient in recipients) {
      
      donor_values <-
        X_scaled[donors, predictors, drop = FALSE]
      
      recipient_values <-
        X_scaled[recipient, predictors]
      
      
      # Variables observed in both recipient and donor
      common <-
        !is.na(donor_values) &
        matrix(
          !is.na(recipient_values),
          nrow = length(donors),
          ncol = length(predictors),
          byrow = TRUE
        )
      
      n_common <- rowSums(common)
      
      
      # Squared differences
      differences <- sweep(
        donor_values,
        2,
        recipient_values,
        "-"
      )^2
      
      differences[!common] <- 0
      
      
      # Sum of available RF weights for each donor
      weight_sum <- rowSums(
        sweep(common, 2, rf_weights, "*")
      )
      
      weighted_sum <- rowSums(
        sweep(differences, 2, rf_weights, "*")
      )
      
      
      distances <- rep(Inf, length(donors))
      
      # Normal RF-weighted distance
      use_rf <- weight_sum > 0
      
      distances[use_rf] <-
        sqrt(
          weighted_sum[use_rf] /
            weight_sum[use_rf]
        )
      
      
      # If all common predictors have RF weight zero,
      # use equal weights across those predictors
      use_equal <- weight_sum == 0 & n_common > 0
      
      distances[use_equal] <-
        sqrt(
          rowSums(
            differences[use_equal, , drop = FALSE]
          ) /
            n_common[use_equal]
        )
      
      
      # Select k nearest valid donors
      valid <- is.finite(distances)
      
      if (!any(valid)) {
        
        X_hat[recipient, target] <-
          mean(X_original[donors, target])
        
      } else {
        
        valid_donors <- donors[valid]
        valid_distances <- distances[valid]
        
        k_use <- min(k, length(valid_donors))
        
        order_k <- order(valid_distances)[seq_len(k_use)]
        
        nearest <- valid_donors[order_k]
        nearest_dist <- valid_distances[order_k]
        nearest_values <- X_original[nearest, target]
        
        # If one or more donors have distance 0,
        # average only those exact matches
        if (any(nearest_dist == 0)) {
          
          X_hat[recipient, target] <-
            mean(nearest_values[nearest_dist == 0])
          
        } else {
          
          donor_weights <- 1 / nearest_dist
          donor_weights <- donor_weights / sum(donor_weights)
          
          X_hat[recipient, target] <-
            sum(donor_weights * nearest_values)
        }
      }
    }
  }
  
  
  list(
    X_hat = X_hat,
    weights = weights
  )
}
