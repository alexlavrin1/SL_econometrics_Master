# ==============================================================================
# Part 1 — Assignment Statistical Learning
# Group 6:
# Alexandre Lavrinenko    xxxxxx
# Ensar Tasgin            646820
# Sanne Maasman           644729
# Sebastiaan van Helden   822236
#
# hybrid_kNN():
#   Imputes missing continuous values using k-nearest neighbours.
#
#   weighted = TRUE:
#     Hybrid kNN. Random-forest variable importance determines how strongly
#     each predictor contributes to the distance calculation.
#
#   weighted = FALSE:
#     Unweighted kNN distance. All predictors receive equal importance.
#
#   In both cases, the final target value is obtained using an
#   inverse-distance weighted average of the selected donors.
# ==============================================================================


hybrid_kNN <- function(X, k, ntree = 500, seed = NULL, weighted = TRUE) {
  # weighted = FALSE: equal weights (standard kNN)
  # Convert input to a matrix and store the number of variables
  X <- as.matrix(X)
  p <- ncol(X)
  
  # If the data has no column names, assign names X1, X2, ...
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(p))
  }
  
  # Keep the original incomplete data unchanged.
  # X_hat will contain the final imputations.
  X_original <- X
  X_hat <- X
  
  # TRUE/FALSE matrix indicating which values were originally missing
  missing <- is.na(X)
  
  # Stores the predictor weights separately for each imputed target variable
  weights <- list()
  
  
  # ============================================================================
  # 1. Scaling for the kNN distance
  # ============================================================================
  
  # Calculate the IQR of each variable.
  # Scaling is needed because variables with a larger numerical scale would
  # otherwise have a larger influence on the Euclidean distance.
  scales <- apply(X, 2, IQR, na.rm = TRUE)
  
  # If a variable has IQR = 0, use its standard deviation instead
  zero_iqr <- scales == 0
  
  scales[zero_iqr] <-
    apply(X[, zero_iqr, drop = FALSE], 2, sd, na.rm = TRUE)
  
  # If a variable is constant, both IQR and SD can be zero.
  # Setting the scale to 1 leaves that variable unchanged.
  scales[is.na(scales) | scales == 0] <- 1
  
  # Divide every column by its scale.
  # No centering is required because subtracting a common mean/median
  # does not affect pairwise differences.
  X_scaled <- sweep(X, 2, scales, "/")
  
  
  # ============================================================================
  # 2. Random-forest variable-importance weights
  # ============================================================================
  
  # For one target variable, fit a separate random forest and return
  # normalized permutation-importance weights for all remaining predictors.
  get_rf_weights <- function(target) {
    
    # All variables except the target are used as predictors
    predictors <- setdiff(seq_len(p), target)
    predictor_names <- colnames(X)[predictors]
    
    # The random forest can only be trained on rows where the target
    # variable itself is observed
    rf_rows <- !is.na(X[, target])
    
    # Response and predictor matrix for the random forest
    y_rf <- X[rf_rows, target]
    x_rf <- X[rf_rows, predictors, drop = FALSE]
    
    
    # --------------------------------------------------------------------------
    # Temporary median imputation for RF predictors
    # --------------------------------------------------------------------------
    
    # Predictors may themselves contain missing values.
    # The randomForest implementation requires complete predictor values,
    # so calculate the median of every predictor.
    medians <- apply(
      X[, predictors, drop = FALSE],
      2,
      median,
      na.rm = TRUE
    )
    
    # Replace missing RF predictor values with their variable median.
    # These temporary values are used ONLY to estimate RF importance;
    # they are not used in the final kNN distance or returned imputations.
    for (j in seq_along(predictors)) {
      x_rf[is.na(x_rf[, j]), j] <- medians[j]
    }
    
    
    # If a seed is supplied, use a different but reproducible seed
    # for each target variable
    if (!is.null(seed)) {
      set.seed(seed + target)
    }
    
    
    # --------------------------------------------------------------------------
    # Fit random forest
    # --------------------------------------------------------------------------
    
    # Fits a regression random forest for:
    # target variable ~ all remaining variables
    fit <- randomForest::randomForest(
      x = x_rf,
      y = y_rf,
      ntree = ntree,
      importance = TRUE
    )
    
    
    # Obtain the default scaled permutation importance (type = 1)
    importance <- drop(
      randomForest::importance(fit, type = 1)
    )
    
    # Attach predictor names to the importance values
    names(importance) <- predictor_names
    
    
    # --------------------------------------------------------------------------
    # Convert RF importance into valid distance weights
    # --------------------------------------------------------------------------
    
    # Negative weights are not meaningful in a distance measure,
    # therefore negative importance values are truncated at zero
    importance <- pmax(importance, 0)
    
    # If no predictor has positive importance, fall back to equal weights
    if (sum(importance) == 0) {
      importance[] <- 1
    }
    
    # Normalize the importance values so that all predictor weights sum to 1
    importance / sum(importance)
  }
  
  
  # ============================================================================
  # 3. kNN imputation
  # ============================================================================
  
  # Repeat the procedure separately for every variable in the dataset
  for (target in seq_len(p)) {
    
    # Recipients are observations where the current target is missing
    recipients <- which(missing[, target])
    
    # If this target has no missing values, move to the next variable
    if (length(recipients) == 0) {
      next
    }
    
    # Donors must have an originally observed value for the target
    donors <- which(!missing[, target])
    
    # Variables used to calculate similarity:
    # all variables except the current target
    predictors <- setdiff(seq_len(p), target)
    
    
    # --------------------------------------------------------------------------
    # Choose predictor weights
    # --------------------------------------------------------------------------
    
    if (weighted) {
      
      # Hybrid kNN:
      # use target-specific random-forest importance weights
      rf_weights <- get_rf_weights(target)
      
    } else {
      
      # Unweighted kNN:
      # all predictors receive equal weight
      rf_weights <- rep(
        1 / length(predictors),
        length(predictors)
      )
      
      names(rf_weights) <- colnames(X)[predictors]
    }
    
    # Store the weights used for this target variable
    weights_name <- colnames(X)[target]
    weights[[weights_name]] <- rf_weights
    
    
    # ==========================================================================
    # 4. Calculate distances and impute each recipient
    # ==========================================================================
    
    for (recipient in recipients) {
      
      # Scaled predictor values for all potential donors
      donor_values <-
        X_scaled[donors, predictors, drop = FALSE]
      
      # Scaled predictor values for the current recipient
      recipient_values <-
        X_scaled[recipient, predictors]
      
      
      # ------------------------------------------------------------------------
      # Determine which predictors are jointly observed
      # ------------------------------------------------------------------------
      
      # A predictor may only contribute to a donor-recipient distance
      # when it is observed for BOTH observations.
      #
      # Each row of 'common' corresponds to one donor and each column
      # corresponds to one predictor.
      common <-
        !is.na(donor_values) &
        matrix(
          !is.na(recipient_values),
          nrow = length(donors),
          ncol = length(predictors),
          byrow = TRUE
        )
      
      # Number of jointly observed predictors for each donor-recipient pair
      n_common <- rowSums(common)
      
      
      # ------------------------------------------------------------------------
      # Calculate squared differences
      # ------------------------------------------------------------------------
      
      # Computes (donor value - recipient value)^2 for every donor
      # and every predictor simultaneously
      differences <- sweep(
        donor_values,
        2,
        recipient_values,
        "-"
      )^2
      
      # Predictors that are not jointly observed should not contribute
      # to the distance
      differences[!common] <- 0
      
      
      # ------------------------------------------------------------------------
      # Renormalize predictor weights for each donor-recipient pair
      # ------------------------------------------------------------------------
      
      # Sum the weights belonging to predictors that are actually available
      # for each donor-recipient comparison.
      #
      # This is required because different donor-recipient pairs may have
      # different sets of jointly observed predictors.
      weight_sum <- rowSums(
        sweep(common, 2, rf_weights, "*")
      )
      
      # Weighted sum of squared differences for every donor
      weighted_sum <- rowSums(
        sweep(differences, 2, rf_weights, "*")
      )
      
      
      # Initialize every donor distance as invalid
      distances <- rep(Inf, length(donors))
      
      
      # ------------------------------------------------------------------------
      # RF-weighted / equal-weight distance
      # ------------------------------------------------------------------------
      
      # Donors for which at least one jointly observed predictor
      # has a positive weight
      use_rf <- weight_sum > 0
      
      # Weighted Euclidean distance.
      # Dividing by weight_sum renormalizes the weights of the predictors
      # actually available for this donor-recipient pair so they sum to 1.
      distances[use_rf] <-
        sqrt(
          weighted_sum[use_rf] /
            weight_sum[use_rf]
        )
      
      
      # ------------------------------------------------------------------------
      # Fallback when all available predictor weights are zero
      # ------------------------------------------------------------------------
      
      # A donor and recipient may share predictors, while all those predictors
      # happen to have RF importance zero.
      use_equal <- weight_sum == 0 & n_common > 0
      
      # In that case, calculate an equal-weight Euclidean distance
      # across the jointly observed predictors
      distances[use_equal] <-
        sqrt(
          rowSums(
            differences[use_equal, , drop = FALSE]
          ) /
            n_common[use_equal]
        )
      
      
      # ==========================================================================
      # 5. Select nearest donors and aggregate their target values
      # ==========================================================================
      
      # Keep only donors for which a valid distance could be calculated
      valid <- is.finite(distances)
      
      if (!any(valid)) {
        
        # If there is no valid donor distance at all,
        # fall back to the observed mean of the target variable
        X_hat[recipient, target] <-
          mean(X_original[donors, target])
        
      } else {
        
        # Extract valid donors and their corresponding distances
        valid_donors <- donors[valid]
        valid_distances <- distances[valid]
        
        # If fewer than k valid donors exist, use all available donors
        k_use <- min(k, length(valid_donors))
        
        # Order donors by increasing distance and select the k nearest
        order_k <- order(valid_distances)[seq_len(k_use)]
        
        nearest <- valid_donors[order_k]
        nearest_dist <- valid_distances[order_k]
        
        # Observed target values belonging to the selected nearest donors
        nearest_values <- X_original[nearest, target]
        
        
        # ----------------------------------------------------------------------
        # Inverse-distance weighted aggregation
        # ----------------------------------------------------------------------
        
        # If one or more selected donors have distance zero, these observations
        # are exact matches according to the available predictors.
        #
        # 1 / distance cannot be calculated for distance = 0, so only
        # the zero-distance donors are averaged in this special case.
        if (any(nearest_dist == 0)) {
          
          X_hat[recipient, target] <-
            mean(nearest_values[nearest_dist == 0])
          
        } else {
          
          # Closer donors receive larger aggregation weights
          donor_weights <- 1 / nearest_dist
          
          # Normalize donor weights so they sum to 1
          donor_weights <- donor_weights / sum(donor_weights)
          
          # Final imputation:
          # inverse-distance weighted average of the k selected donor values
          X_hat[recipient, target] <-
            sum(donor_weights * nearest_values)
        }
      }
    }
  }
  
  
  # ============================================================================
  # 6. Output
  # ============================================================================
  
  # X_hat:
  #   original data with missing values replaced by kNN imputations
  #
  # weights:
  #   predictor weights used separately for every incomplete target variable
  list(
    X_hat = X_hat,
    weights = weights
  )
}
