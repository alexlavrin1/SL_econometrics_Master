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
library(ggplot2)
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
# strength: controls how strongly missingness probability depends on the driver in MAR/MNAR
make_missing <- function(X, vars, prop = 0.2,
                         mechanism = c("MCAR", "MAR", "MNAR"),
                         driver = NULL, strength = 2) {
  
  # Validates that mechanism is one of the allowed options; defaults to MCAR if not specified
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
                    MCAR = rep(0, n),
                    MAR  = strength * z_driver,
                    MNAR = strength * z_driver + strength * as.numeric(scale(X[, j]))
    )
    
  
    # Draws m row indices out of n (without replacement): probability of each being selected is (roughly) proportional to plogis(score) for that row
    idx <- sample(n, m, prob = plogis(score))
    
    # Sets value NA to those m rows in column j
    X_miss[idx, j] <- NA
  }
  
  # Returns the full matrix with missing values (NA) now inserted across the vars columns
  X_miss
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

#addition: bias. NRMSE tells how mich it is off, bias says in which direction
#report per variable which is imputed 

evaluate <- function(X_true, X_imp, X_miss) {
  miss_mask <- is.na(X_miss)                 # TRUE/FALSE matrix (same size as X_miss)
  diff      <- X_imp - X_true                # elementwise error matrix; 0 on non-missing (observed) cells 
  vars      <- which(colSums(miss_mask) > 0) # collects the indices of columns with missing values
  
  # Assigns column names if no column names
  var_names <- colnames(X_true)
  if (is.null(var_names)) var_names <- paste0("X", seq_len(ncol(X_true)))
  
  # Applies function(j) once to each value of j in vars
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

source("hybrid_kNN.R")
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
saveRDS(list(results = results, weights = weights_df,
             settings = list(n = n, p = p, k = k, prop = prop, strength = strength,
                             miss_vars = miss_vars, driver = driver, R = R,
                             base_seed = base_seed)),
        "sim_results.rds")


#=====================================================================================================================================
#=====================================================================================================================================
# 2.5 Summarize and visualize data data  ------------------------------------------------------
#plot before and after imputations
#=====================================================================================================================================
#=====================================================================================================================================
## 2.5a Summary statistics ----

mech_levels <- c("MCAR", "MAR", "MNAR")

# 1. Mean NRMSE and bias per scenario x mechanism x method (averaged over reps and X1-X3)
tab <- aggregate(cbind(nrmse, bias) ~ scenario + mechanism + method,
                 data = results, FUN = mean)
tab$mechanism <- factor(tab$mechanism, levels = mech_levels)
tab <- tab[order(tab$scenario, tab$mechanism, tab$method), ]

cat("\n===== Mean NRMSE (lower = better) =====\n")
print(ftable(round(xtabs(nrmse ~ scenario + mechanism + method, data = tab), 3)))

cat("\n===== Mean bias (0 = unbiased) =====\n")
print(ftable(round(xtabs(bias ~ scenario + mechanism + method, data = tab), 3)))


# 2. Paired difference hybrid - standard (per rep, averaged over X1-X3 first)
per_rep <- aggregate(nrmse ~ scenario + mechanism + method + rep,
                     data = results, FUN = mean)

d <- merge(subset(per_rep, method == "hybrid"),
           subset(per_rep, method == "unweighted"),
           by = c("scenario", "mechanism", "rep"),
           suffixes = c("_hyb", "_std"))
d$diff <- d$nrmse_hyb - d$nrmse_std

paired <- do.call(rbind, lapply(split(d, list(d$scenario, d$mechanism), drop = TRUE),
                                function(x) {
                                  se <- sd(x$diff) / sqrt(nrow(x))
                                  data.frame(scenario      = x$scenario[1],
                                             mechanism     = x$mechanism[1],
                                             mean_diff     = mean(x$diff),
                                             se            = se,
                                             t_stat        = mean(x$diff) / se,
                                             hybrid_better = mean(x$diff < 0))   # share of reps where hybrid wins
                                }))
paired$mechanism <- factor(paired$mechanism, levels = mech_levels)
paired <- paired[order(paired$scenario, paired$mechanism), ]
rownames(paired) <- NULL

cat("\n===== Paired difference NRMSE (hybrid - standard; < 0 = hybrid better) =====\n")
print(transform(paired,
                mean_diff = round(mean_diff, 4),
                se        = round(se, 4),
                t_stat    = round(t_stat, 2)))


# 3. Average RF weights when imputing X1 (pooled over mechanisms)
w_avg <- aggregate(weight ~ scenario + predictor,
                   data = subset(weights_df, target == "X1"), FUN = mean)
w_avg$predictor <- factor(w_avg$predictor, levels = paste0("X", 2:p))

cat("\n===== Average RF weights for imputing X1 =====\n")
print(round(xtabs(weight ~ scenario + predictor, data = w_avg), 3))


# 4. Save tables for the report
write.csv(tab,    "table_nrmse_bias.csv",   row.names = FALSE)
write.csv(paired, "table_paired_diff.csv",  row.names = FALSE)


install.packages("VIM", dependencies = TRUE)
library(VIM)

illus_rep <- 1   # same replication as rep 1 in the simulation

## One VIM marginplot: X4 (driver) vs X1 (imputed), saved as pdf ----
vim_plot <- function(scenario, mechanism, method, file) {
  set.seed(base_seed + illus_rep)          # same data and same missing pattern for every method
  X      <- gen_data(scenario)
  X_miss <- make_missing(X, miss_vars, prop, mechanism,
                         driver = driver, strength = strength)
  X_imp  <- methods[[method]](X_miss)$X_hat
  
  # VIM needs indicator columns "<var>_imp" that mark which values were imputed
  df <- data.frame(
    X4     = X_imp[, "X4"],
    X1     = X_imp[, "X1"],
    X4_imp = is.na(X_miss[, "X4"]),        # always FALSE (driver is fully observed)
    X1_imp = is.na(X_miss[, "X1"])
  )
  
  pdf(file, width = 4, height = 4)
  marginplot(df, delimiter = "_imp")
  dev.off()
}

## 2 scenarios x 3 mechanisms x 2 methods = 12 pdfs ----
for (s in c("A", "B")) {
  for (mech in c("MCAR", "MAR", "MNAR")) {
    for (m in c("unweighted", "hybrid")) {
      vim_plot(s, mech, m, sprintf("vim_%s_%s_%s.pdf", s, mech, m))
    }
  }
}

## NRMSE per situation for the subcaptions (average over R replications, X1) ----
aggregate(nrmse ~ scenario + mechanism + method,
          data = subset(results, variable == "X1"), FUN = mean)



#to see true vs imputed values

illus_rep <- 1

## True vs imputed values of X1 for one replication ----
make_tvi <- function(scenario, mechanism) {
  set.seed(base_seed + illus_rep)
  X      <- gen_data(scenario)
  X_miss <- make_missing(X, miss_vars, prop, mechanism,
                         driver = driver, strength = strength)
  miss   <- is.na(X_miss[, "X1"])
  
  out <- lapply(c("hybrid", "unweighted", "mean"), function(m) {
    X_imp <- methods[[m]](X_miss)$X_hat
    data.frame(true = X[miss, "X1"], imputed = X_imp[miss, "X1"], method = m)
  })
  df <- do.call(rbind, out)
  df$scenario  <- scenario
  df$mechanism <- mechanism
  df
}

combos <- expand.grid(scenario  = c("A", "B"),
                      mechanism = c("MCAR", "MAR", "MNAR"),
                      stringsAsFactors = FALSE)

tvi <- do.call(rbind, lapply(seq_len(nrow(combos)), function(i) {
  make_tvi(combos$scenario[i], combos$mechanism[i])
}))

tvi$mechanism <- factor(tvi$mechanism, levels = c("MCAR", "MAR", "MNAR"))
tvi$scenario  <- factor(tvi$scenario, levels = c("A", "B"),
                        labels = c("Scenario A", "Scenario B"))
tvi$method    <- factor(tvi$method, levels = c("hybrid", "unweighted", "mean"),
                        labels = c("Hybrid kNN", "Standard kNN", "Mean"))

p_tvi <- ggplot(tvi, aes(x = true, y = imputed, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.5, size = 1) +
  facet_grid(scenario ~ mechanism) +
  coord_equal() +
  scale_colour_manual(values = c("Hybrid kNN"   = "#E69F00",
                                 "Standard kNN" = "#CC3311",
                                 "Mean"         = "grey60")) +
  labs(x = expression("True value of " * X[1]),
       y = expression("Imputed value of " * X[1]),
       colour = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_tvi)
ggsave("fig_true_vs_imputed.pdf", p_tvi, width = 9, height = 6)
