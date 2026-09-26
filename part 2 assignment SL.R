# ==============================================================================
# Part 2 — Assignment Statistical Learning
# Group 6:
# Alexandre Lavrinenko    xxxxxx
# Ensar Tasgin:           646820
# Sanne Maasman           xxxxxx 
# Sebastiaan van Helden   XXXXX

#Build up as follows:
#2.1: Builds data generating process
#2.2: builds the missing data mechanisms

#remaining part run the simulation and comparison
#2.3: makes scoring fucntion to evaluate the imputations
#2.4: Runs the simulation
#2.5: Summarizes and plots the results
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
# A: equal correlation --> normal kNN should perform better --> see report 
# B: heabily correlated variables --> hybrid should perform better --> see report
  # with 3 variables heavily correlated 

generate_data_A <- function(n, p = 10, rho = 0.5) {
  Sigma <- matrix(rho, nrow = p, ncol = p)
  diag(Sigma) <- 1
  
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("X", seq_len(p))
  X
}


generate_data_B <- function(n, p = 10, n_signal = 4, var_e = 0.3) {
  Z <- rnorm(n, mean = 0, sd = 1)
  
  signal <- replicate(n_signal, Z + rnorm(n, sd = sqrt(var_e)))
  
  noise <- matrix(rnorm(n * (p - n_signal), mean = 0, sd = 1),
                  nrow = n, ncol = p - n_signal)
  
  X <- cbind(signal, noise)
  colnames(X) <- paste0("X", seq_len(p))
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
## Probabilities via a logistic function on standardized values, so the
## mechanism is scale-independent and also works on unseen data.
## sample(n, m, prob = w) guarantees exactly m missing values per variable.

make_missing <- function(X, vars, prop = 0.2,
                         mechanism = c("MCAR", "MAR", "MNAR"),
                         driver = NULL, strength = 2) {
  
  mechanism <- match.arg(mechanism)
  n <- nrow(X)
  m <- round(prop * n)
  X_miss <- X
  
  if (mechanism %in% c("MAR", "MNAR")) {
    if (is.null(driver)) {
      stop("driver just be specified MAR/MNAR")
    }
    if (driver %in% vars) {
      stop("driver must be fully observed: choose a column besued 'vars'")
    }
    z_driver <- as.numeric(scale(X[, driver]))
  }
  
  for (j in vars) {
    score <- switch(mechanism,
                    MCAR = rep(0, n),
                    MAR  = strength * z_driver,
                    MNAR = strength * z_driver + strength * as.numeric(scale(X[, j]))
    )
    
    idx <- sample(n, m, prob = plogis(score))
    X_miss[idx, j] <- NA
  }
  
  X_miss
}


## Check: does the mechanism actually work ----
##   MCAR : all means approximately equal
##   MAR  : driver mean differs between missing and observed rows
##   MNAR : additionally, the mean of the variable itself differs

check_missing <- function(X, X_miss, vars, driver = NULL) {
  for (j in vars) {
    r <- is.na(X_miss[, j])
    cat(sprintf("Var %-3d | prop = %.3f", j, mean(r)))
    if (!is.null(driver)) {
      cat(sprintf(" | driver mis/obs = %6.2f / %6.2f",
                  mean(X[r, driver]), mean(X[!r, driver])))
    }
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
#     1 measn equally as good as imputing average under MCAR
# assignment asks which evaluation metrix do we suggest --> choose 1: NRMSE

#maybe addition: bias. NRMSE tells how mich it is off, bias says in which direction
#rapport per variable which is imputed 

evaluate <- function(X_true, X_imp, X_miss) {
  miss_mask <- is.na(X_miss)                 # where we made NA's
  diff      <- X_imp - X_true                # elementwise; 0 on ibserved cells
  vars      <- which(colSums(miss_mask) > 0) # only columns with missing values
  
  var_names <- colnames(X_true)
  if (is.null(var_names)) var_names <- paste0("X", seq_len(ncol(X_true)))
  
  nrmse <- vapply(vars, function(j) {
    e <- diff[miss_mask[, j], j]             # only missing cells from column j
    sqrt(mean(e^2)) / sd(X_true[, j])
  }, numeric(1))
  
  bias <- vapply(vars, function(j) {
    mean(diff[miss_mask[, j], j])
  }, numeric(1))
  
  data.frame(variable = var_names[vars], nrmse = nrmse, bias = bias,
             row.names = NULL)
}

## Baseline + sanity check: average imputation ----
##  Expectation under MCAR: NRMSE ~ 1, bias ~ 0
mean_impute <- function(X_miss) {
  X_imp <- X_miss
  for (j in seq_len(ncol(X_miss))) {
    na <- is.na(X_miss[, j])
    X_imp[na, j] <- mean(X_miss[, j], na.rm = TRUE)
  }
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
#   - A: regular kNN should perform better
#   - B: hybrid kNN should perform better
# 
# Then per set, 4 matrices will be made:
#   1. original (complete) matrix
#   2. matirx with deleted values by MCAR, MAR & MNAR
#   3. matrix with imputed valued based on matrix 2
#   4. differences matrix where matrix with imputed values - orginal matrix (to evaluate imputations, omly imputations have nonzero values (besides perfect imputations))
#   
# Them evaluation function is used on differneces matrix
# Results will be plotted (2.5)

source("hybrid_kNN v2 (not OG).R")
args(hybrid_kNN)   

## Settings (motivate them in report)
n         <- 500      #num obs
p         <- 10       #num vars
k         <- 5        # num donors for kNN
prop      <- 0.2      #proportion missing
strength  <- 2        #how mucht the missingness depends on the dpeendent variable in MNAR, MAR
miss_vars <- 1:3      # in which columns missing vairables will be present could also be: c(1, 2, 7) om ook een ruisvariabele te imputeren in B
driver    <- 4        # The variable which decides which obs are going to be missing (MNAR, MAR)
R         <- 100      # number of repitiions per scenario and mechanism (needed for noise in randomness of variables used) (see if we keep this in)
base_seed <- 2026     #seed 


## Data per scenario ----
gen_data <- function(scenario) {
  switch(scenario,
         A = generate_data_A(n, p),
         B = generate_data_B(n, p))
}


methods <- list(
  hybrid     = function(Xm) hybrid_kNN(Xm, k = k, weighted = TRUE),
  unweighted = function(Xm) hybrid_kNN(Xm, k = k, weighted = FALSE),
  mean       = function(Xm) list(X_hat = mean_impute(Xm), weights = NULL)
)

## RF-weight to table for plotting in 2.5 ----
weights_to_df <- function(w) {
  do.call(rbind, lapply(names(w), function(t) {
    data.frame(target = t, predictor = names(w[[t]]),
               weight = unname(w[[t]]))
  }))
}

## one repetition
one_run <- function(scenario, mechanism, rep) {
  
  set.seed(base_seed + rep)   #same seed per repetition
  
  X      <- gen_data(scenario)
  X_miss <- make_missing(X, miss_vars, prop, mechanism,
                         driver = driver, strength = strength)
  
  scores  <- list()
  weights <- NULL
  
  for (m in names(methods)) {
    out   <- methods[[m]](X_miss)            # every method on the same X_miss
    ev    <- evaluate(X, out$X_hat, X_miss)
    ev$method   <- m
    scores[[m]] <- ev
    
    if (m == "hybrid") weights <- weights_to_df(out$weights)
  }
  
  scores <- do.call(rbind, scores)
  scores$scenario  <- scenario
  scores$mechanism <- mechanism
  scores$rep       <- rep
  
  weights$scenario  <- scenario
  weights$mechanism <- mechanism
  weights$rep       <- rep
  
  list(scores = scores, weights = weights)
}

## Grid: 2 scenarios x 3 mechanisms x R repitions ----
grid <- expand.grid(scenario  = c("A", "B"),
                    mechanism = c("MCAR", "MAR", "MNAR"),
                    rep       = seq_len(R),
                    stringsAsFactors = FALSE)

## running
t0 <- Sys.time()
runs <- lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  cat(sprintf("%s | %s | %-4s | rep %d\n",
              format(Sys.time(), "%H:%M:%S"), g$scenario, g$mechanism, g$rep))
  one_run(g$scenario, g$mechanism, g$rep)
})
print(Sys.time() - t0)

results    <- do.call(rbind, lapply(runs, `[[`, "scores"))
weights_df <- do.call(rbind, lapply(runs, `[[`, "weights"))

saveRDS(list(results = results, weights = weights_df,
             settings = list(n = n, p = p, k = k, prop = prop, strength = strength,
                             miss_vars = miss_vars, driver = driver, R = R,
                             base_seed = base_seed)),
        "sim_results.rds")


## Sanity checks (set run_checks TRUE to run) ----
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
#do evaluation








#=====================================================================================================================================
#=====================================================================================================================================
# 2.5 Summarize and visualize data data  ------------------------------------------------------
#plot before and after imputations
#=====================================================================================================================================
#=====================================================================================================================================








