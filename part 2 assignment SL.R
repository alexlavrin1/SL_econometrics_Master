# ==============================================================================
# Part 2 — Assignment Statistical Learning
# Group 6:
# Alexandre Lavrinenko    xxxxxx
# Ensar Tasgin:           646820
# Sanne Maasman           644729 
# Sebastiaan van Helden   XXXXX

# Build up as follows:
# 2.1: Builds data generating process
# 2.2: Builds the missing data mechanisms

# Remaining part run the simulation and comparison
# 2.3: Makes scoring function (evaluation metric) to evaluate the imputations
# 2.4: Runs the simulation
# 2.5: Summarizes and plots the results
# ==============================================================================

# Libraries ---------------------------------------------------------------

library(MASS)
set.seed(2026)

#=====================================================================================================================================
#=====================================================================================================================================
# 2.1 Data generating processes -------------------------------------------
#=====================================================================================================================================
#=====================================================================================================================================
# Two scenarios:
# A: equal correlation --> normal kNN should perform better (should they not perform equally well?) --> see report 
# B: heavily correlated variables --> hybrid should perform better --> see report
  # with 3 variables heavily correlated 

generate_data_A <- function(n, p = 10, rho = 0.5) {
  # Builds a pxp matrix filled entirely with rho (0.5)
  Sigma <- matrix(rho, nrow = p, ncol = p)
  
  # Overwrites the diagonal with 1's
  diag(Sigma) <- 1
  
  # Draws n observations from a p-dimensional multivariate normal with mean 0 for every variable and correlation matrix Sigma
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  # Names columns X1,...,X10
  colnames(X) <- paste0("X", seq_len(p))
  
  # Returns an nxp matrix
  X
}

generate_data_B <- function(n, p = 10, n_signal = 4, var_e = 0.3) {
  # Generates a latent common factor (standard normal, vector length n)
  Z <- rnorm(n, mean = 0, sd = 1)
  
  # Generates an nxn_signal matrix, with each column Z plus an independent noise term with variance var_e 
  signal <- replicate(n_signal, Z + rnorm(n, sd = sqrt(var_e)))
  
  # Draws n * (p - n_signal) random values from a standard normal distribution, then reshapes these values into a nx(p-n_signal) matrix
  noise <- matrix(rnorm(n * (p - n_signal), mean = 0, sd = 1), 
                  nrow = n, ncol = p - n_signal)

  # Combines matrices signal and noise into one large matrix
  X <- cbind(signal, noise)

  # Assigns column names X1…X10
  colnames(X) <- paste0("X", seq_len(p))

  # Returns an nxp matrix
  X
}



#=====================================================================================================================================
#=====================================================================================================================================
# 2.2 Missing data mechanisms ---------------------------------------------
#=====================================================================================================================================
#=====================================================================================================================================
## MCAR:  P(R | X) = P(R)
## MAR:   P(R | X) = P(R | X_obs)
## MNAR:  P(R | X) = P(R | X_obs, X_mis)
##
## Probabilities via a logistic function on standardized values, so the mechanism is scale-independent and also works on unseen data.
## sample(n, m, prob = w) guarantees exactly m missing values per variable.


# vars: specifies in which columns we induce missing values
# prop: target proportion missing per column
# driver: the fully-observed variable that missingness probability depends on under MAR/MNAR
# strength: controls how strongly missingness probability depends on the driver in MAR/MNAR (and for MNAR, on the variable's own value too)
make_missing <- function(X, vars, prop = 0.2,
                         mechanism = c("MCAR", "MAR", "MNAR"),
                         driver = NULL, strength = 2) {
  
  # Validates that mechanism is one of the allowed options; defaults to MCAR if not specified and throws error if user chose a non-existing parameter input
  mechanism <- match.arg(mechanism)
  
  # Number of rows
  n <- nrow(X)
  
  # Number of cells to delete per column, rounded to nearest integer (for example, prop=0.2 and n=500 imply m=100).
  m <- round(prop * n)
  
  # Copy of X, such that the original X stays unchanged
  X_miss <- X

  
  if (mechanism %in% c("MAR", "MNAR")) {
    if (is.null(driver)) {
      stop("driver just be specified MAR/MNAR")
    }
    if (driver %in% vars) {
      stop("driver must be fully observed: choose a column besued 'vars'")
    }
    # Standardizes the driver column such that 'strength' has a comparable effect regardless of the driver's original scale
    z_driver <- as.numeric(scale(X[, driver]))
  }

  # Iterates over every column in 'vars'
  for (j in vars) {
    # For each observation, assigns a score; rows with higher scores are more likely to be chosen as missing.
    score <- switch(mechanism,
                    # Identical probability of missing for all observations (score = 0 for each) 
                    MCAR = rep(0, n),
                    # Score depends linearly on the (standardized) driver only
                    MAR  = strength * z_driver,
                    # Score depends on the (standardized) driver and the standardized value of the variable being deleted; its own value influences whether it goes missing.
                    MNAR = strength * z_driver + strength * as.numeric(scale(X[, j]))
    )

    
    # plogis(score): converts the score to a probability in (0,1) using the logistic (sigmoid) function
    
    # Draws m row indices out of n (without replacement): probability of each being selected is (roughly) proportional to plogis(score) for that row
    idx <- sample(n, m, prob = plogis(score))

    # Sets value NA to those m rows in column j
    X_miss[idx, j] <- NA
  }

  # Returns the full matrix with missing values (NA) now inserted across the vars columns
  X_miss
}




## Diagnostic function: confirm whether the mechanisms actually work
##   MCAR : all means approximately equal
##   MAR  : driver mean differs between missing and observed rows
##   MNAR : additionally, the mean of the variable itself differs

check_missing <- function(X, X_miss, vars, driver = NULL) {
  for (j in vars) {
    # Returns a vector of TRUE/FALSE values, one for each row in column j (each entry tells you whether that specific cell is missing (NA) or not)
    r <- is.na(X_miss[, j])

    # Prints Var j | prop = mean(r), where mean(r) is the actual observed proportion of missing values in that column
    cat(sprintf("Var %-3d | prop = %.3f", j, mean(r)))

    # If a driver was specified, print:
          # The average driver value among observations where variable j went missing,
          # The average driver value among observations where variable j is still observed.
    if (!is.null(driver)) {
      cat(sprintf(" | driver mis/obs = %6.2f / %6.2f",
                  mean(X[r, driver]), mean(X[!r, driver])))
    }

    # X[r,j]: from column j of the true data X, give me only the rows where r is TRUE" (that is, only the rows where variable j was deleted in X_miss). 
    # This recovers what those now-missing values actually were, before they got wiped out.
    # mean(X[r,j]): averages these values
    
    # Prints:
          # The true average of variable j, restricted to rows that ended up missing,
          # The true average of variable j, restricted to rows that stayed observed.
    cat(sprintf(" | self mis/obs = %6.2f / %6.2f\n",
                mean(X[r, j]), mean(X[!r, j])))
  }
}







#=====================================================================================================================================
#=====================================================================================================================================
# 2.3 Evaluation metric ---------------------------------------------------
#=====================================================================================================================================
#=====================================================================================================================================
# NRMSE: Normalized Root Mean Squared Error (RMSE/stdev())
#     0 means perfect imputation
#     1 means equally as good as imputing average under MCAR
# assignment asks which evaluation metric do we suggest --> choose 1: NRMSE

#maybe addition: bias. NRMSE tells how mich it is off, bias says in which direction
#report per variable which is imputed 

evaluate <- function(X_true, X_imp, X_miss) {
  miss_mask <- is.na(X_miss)                 # TRUE/FALSE matrix (same size as X_miss)
  diff      <- X_imp - X_true                # elementwise error matrix; 0 on non-missing (observed) cells 
  vars      <- which(colSums(miss_mask) > 0) # collects the indices of columns with missing values

  # Assigns column names if no column names
  var_names <- colnames(X_true)
  if (is.null(var_names)) var_names <- paste0("X", seq_len(ncol(X_true)))

  # Applies function(j) once to each value of j in vars (j=1, then j=2, then j=3...)
  nrmse <- vapply(vars, function(j) {
                            # Computes imputation errors only for the missing cells from column j
                            e <- diff[miss_mask[, j], j]             
                            sqrt(mean(e^2)) / sd(X_true[, j])
                        # Returns a single number
                        }, numeric(1))

  # Computes bias only for the rows in column j where data was missing, then averages those errors (mean signed error)
  bias <- vapply(vars, function(j) {
    mean(diff[miss_mask[, j], j])
  }, numeric(1))

  # Collects all computed values into a dataframe
  data.frame(variable = var_names[vars], nrmse = nrmse, bias = bias,
             row.names = NULL)
}



## Baseline + sanity check: average imputation
##  Expectation under MCAR: NRMSE ~ 1, bias ~ 0
mean_impute <- function(X_miss) {
  X_imp <- X_miss
  # Iterate over every column in X_miss
  for (j in seq_len(ncol(X_miss))) {
    # Returns which rows are missing in column j (TRUE/FALSE per row)
    na <- is.na(X_miss[, j])
    # Fill missing rows in column j with that column's own observed mean (na.rm = TRUE excludes the NAs themselves when computing the mean)
    X_imp[na, j] <- mean(X_miss[, j], na.rm = TRUE)
  }
  # Returns the fully imputed matrix 
  X_imp
}




#=====================================================================================================================================
#=====================================================================================================================================
# 2.4 Simulation kNN ------------------------------------------------------
#=====================================================================================================================================
#=====================================================================================================================================

# How is the simulatiom going to go?
#   
# First, two datasets are generates:
#   - A: regular kNN should perform well.
#   - B: hybrid kNN should perform well.
# 
# Then per set, 4 matrices will be made:
#   1. original (complete) matrix
#   2. matrix with deleted values by MCAR, MAR & MNAR
#   3. matrix with imputed valued based on matrix 2
#   4. differences matrix is constructed where matrix with imputed values - orginal matrix (to evaluate imputations, only imputations have nonzero values (besides perfect imputations))
#   
# Then evaluation function is used on differences matrix
# Results will be plotted in Section (2.5)

source("hybrid_kNN v2 (not OG).R")
args(hybrid_kNN)   

## Settings (motivate them in report)
n         <- 500      # num obs
p         <- 10       # num vars
k         <- 5        # num donors for kNN
prop      <- 0.2      # proportion missing
strength  <- 2        # how much the missingness depends on the dependent variable in MNAR, MAR
miss_vars <- 1:3      # in which columns missing variables will be present could also be: c(1, 2, 7) om ook een ruisvariabele te imputeren in B
driver    <- 4        # The variable which decides which obs are going to be missing (MNAR, MAR)
R         <- 100      # Number of repitiions per scenario and mechanism (needed for noise in randomness of variables used) (see if we keep this in)
base_seed <- 2026     #seed 


## Generating the data for each scenario
gen_data <- function(scenario) {
  switch(scenario,
         A = generate_data_A(n, p),
         B = generate_data_B(n, p))
}

# Named list of the three imputation methods being compared, each wrapped as a function of just Xm (the data with missing values)
methods <- list(
  hybrid     = function(Xm) hybrid_kNN(Xm, k = k, weighted = TRUE),
  unweighted = function(Xm) hybrid_kNN(Xm, k = k, weighted = FALSE),
  mean       = function(Xm) list(X_hat = mean_impute(Xm), weights = NULL)
)

## Converts RF-weight to table for easier plotting in 2.5: combines smaller data frames for each target variable (X1, X2, X3) into a larger dataframe
weights_to_df <- function(w) {
  do.call(rbind, lapply(names(w), function(t) {
    data.frame(target = t, predictor = names(w[[t]]),
               weight = unname(w[[t]]))
  }))
}

# Defines a function that runs one full repetition of the simulation, given a scenario, missingness mechanism and a repetition number.
one_run <- function(scenario, mechanism, rep) {
  
  # Ensures the same seed per repetition (same seed is reused across different mechanisms for the same rep)
  set.seed(base_seed + rep)   
  
  X      <- gen_data(scenario)
  X_miss <- make_missing(X, miss_vars, prop, mechanism,
                         driver = driver, strength = strength)

  # Intitializes variables 
  scores  <- list()
  weights <- NULL
  
  for (m in names(methods)) {
    out   <- methods[[m]](X_miss)            # every method on the same X_miss
    ev    <- evaluate(X, out$X_hat, X_miss)  # computes NRMSE and bias per affected variable: compares true data X against this method's imputed result out$X_hat
    ev$method   <- m                         # adds a new column (method) to the evaluation data frame that labels every row with which method produced it
    scores[[m]] <- ev                        # stores this method's labelled evaluation data frame into the scores list, under the name m
    
    if (m == "hybrid") weights <- weights_to_df(out$weights)
  }

  scores <- do.call(rbind, scores)            # Combines the scores data frames row-wise into one large dataframe
  scores$scenario  <- scenario                # Adds "scenario" column to scores dataframe
  scores$mechanism <- mechanism               # Adds "mechanism" column to scores dataframe
  scores$rep       <- rep                     # Adds "rep" column to scores dataframe
  
  weights$scenario  <- scenario               # Adds "scenario" column to weights dataframe
  weights$mechanism <- mechanism              # Adds "mechanism" column to weights dataframe
  weights$rep       <- rep                    # Adds "rep" column to weights dataframe
  
  list(scores = scores, weights = weights)    # Returns the scores and weights dataframes (as a list)
}


## Grid: 2 scenarios x 3 mechanisms x R repitions 

# Takes multiple vectors and builds a data frame containing every possible combination of their elements (one row per combination)
grid <- expand.grid(scenario  = c("A", "B"),
                    mechanism = c("MCAR", "MAR", "MNAR"),
                    rep       = seq_len(R),
                    stringsAsFactors = FALSE)

## running
t0 <- Sys.time()
runs <- lapply(seq_len(nrow(grid)), function(i) {       # runs function once for every index in nrow(grid)
  g <- grid[i, ]                                        # extracts row i of grid (scenario, mechanism, rep)
  cat(sprintf("%s | %s | %-4s | rep %d\n",
              format(Sys.time(), "%H:%M:%S"), g$scenario, g$mechanism, g$rep))
  one_run(g$scenario, g$mechanism, g$rep)               # calls one_run() with this row's scenario, mechanism and rep
})
print(Sys.time() - t0)

results    <- do.call(rbind, lapply(runs, `[[`, "scores"))   # Extracts 'scores' from every iteration in runs (one data frame per iteration) and stacks them row-wise into one large dataframe
weights_df <- do.call(rbind, lapply(runs, `[[`, "weights"))  # Same idea, but extracting and stacking the 'weights' element instead

# Saves results, weights, and settings to disk (such that simulation does not need to be rerun every time want to explore/plot results).
saveRDS(list(results = results, weights = weights_df
             settings = list(n = n, p = p, k = k, prop = prop, strength = strength,
                             miss_vars = miss_vars, driver = driver, R = R,
                             base_seed = base_seed)),
        "sim_results.rds")


## Sanity checks (set run_checks TRUE to run)
run_checks <- FALSE
if (run_checks) {
  X_B    <- generate_data_B(n = 2000, p = 10)
  print(round(cor(X_B), 2))
  
  X_mcar <- make_missing(X_B, miss_vars, 0.2, "MCAR")
  X_mar  <- make_missing(X_B, miss_vars, 0.2, "MAR",  driver = driver)
  X_mnar <- make_missing(X_B, miss_vars, 0.2, "MNAR", driver = driver)
  cat("\n--- MCAR ---\n"); check_missing(X_B, X_mcar, miss_vars, driver)
  cat("\n--- MAR  ---\n"); check_missing(X_B, X_mar,  miss_vars, driver)
  cat("\n--- MNAR ---\n"); check_missing(X_B, X_mnar, miss_vars, driver)
  
  # evaluate: mean-imputatie with MCAR must be NRMSE ~ 1 and bias ~ 0 
  print(evaluate(X_B, mean_impute(X_mcar), X_mcar))
  
  # weighted = FALSE must give equal weights
  res_uw <- hybrid_kNN(X_mcar[1:500, ], k = 5, weighted = FALSE)
  print(res_uw$weights$X1)
}



#===============================================
#do evaluation (already done if I'm not mistaken)








#=====================================================================================================================================
#=====================================================================================================================================
# 2.5 Summarize and visualize data data  ------------------------------------------------------
#plot before and after imputations
#=====================================================================================================================================
#=====================================================================================================================================








