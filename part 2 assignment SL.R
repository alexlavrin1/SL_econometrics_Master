# ==============================================================================
# Part 2 — Assignment Statistical Learning
# Group 6:
# Alexandre Lavrinenko    xxxxxx
# Ensar Tasgin:           646820
# Sanne Maasman           xxxxxx 
# Sebastiaan van Helden   822236

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


#to be deleted
## Scenario A: equicorrelated — alle variabelen even informatief ----
##   Alle paarsgewijze correlaties gelijk aan rho.
##   Verwachting: RF-importances vlak -> gewogen ~ ongewogen afstand,
##   hybrid heeft geen voordeel (alleen extra schattingsruis).

generate_data_A <- function(n, p = 10, rho = 0.5) {
  Sigma <- matrix(rho, nrow = p, ncol = p)
  diag(Sigma) <- 1
  
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("X", seq_len(p))
  X
}

#to be deleted
## Scenario B: signaalblok + ruis ----
##   X1..X_signal = Z + e_i  (gedeelde latente factor)
##   X_(signal+1)..Xp        = onafhankelijke ruis
##   Impliciete correlatie binnen het blok: Var(Z) / (Var(Z) + var_e)
##   Verwachting: RF legt bijna al het gewicht op het signaalblok,
##   ongewogen afstand wordt verstoord door de ruisdimensies.

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
## Kansen via logistische functie op gestandaardiseerde waarden, zodat het
## mechanisme schaal-onafhankelijk is en ook werkt op ongeziene data.
## sample(n, m, prob = w) garandeert exact m missings per variabele.

make_missing <- function(X, vars, prop = 0.2,
                         mechanism = c("MCAR", "MAR", "MNAR"),
                         driver = NULL, strength = 2) {
  
  mechanism <- match.arg(mechanism)
  n <- nrow(X)
  m <- round(prop * n)
  X_miss <- X
  
  if (mechanism %in% c("MAR", "MNAR")) {
    if (is.null(driver)) {
      stop("driver moet gespecificeerd zijn voor MAR/MNAR")
    }
    if (driver %in% vars) {
      stop("driver moet volledig geobserveerd zijn: kies een kolom buiten 'vars'")
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


## Controle: bijt het mechanisme daadwerkelijk? ----
##   MCAR : alle gemiddelden ongeveer gelijk
##   MAR  : driver-gemiddelde wijkt af tussen missing en observed
##   MNAR : daarnaast wijkt ook het gemiddelde van de variabele zelf af

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
# Wacht op part 1: aggregatiemethode bepaalt metric (mean -> RMSE, median -> MAE)
# NRMSE: Normalized Root Mean Squared Error (RMSE/stdev())
#     0 means perfect imputation
#     1 measn equally as good as imputing average under MCAR
# assignment asks which evaluation metrix do we suggest --> choose 1: NRMSE

#maybe addition: bias. NRMSE tells how mich it is off, bias says in which direction
#rapport per variable which is imputed 

# Metric: NRMSE = RMSE / sd(ware variabele), alleen over de weggegooide cellen.
#   Part 1 aggregeert met een (1/d-gewogen) gemiddelde -> minimaliseert kwadratische fout -> RMSE
#   0 = perfecte imputatie, ~1 = even goed als gemiddelde-imputatie (onder MCAR)
# Bias = gemiddelde van (geimputeerd - waar): richting van de fout, belangrijk onder MNAR
# Alles per geimputeerde variabele.

evaluate <- function(X_true, X_imp, X_miss) {
  miss_mask <- is.na(X_miss)                 # TRUE waar we iets weggegooid hebben
  diff      <- X_imp - X_true                # elementwise; 0 op geobserveerde cellen
  vars      <- which(colSums(miss_mask) > 0) # alleen kolommen met missings
  
  var_names <- colnames(X_true)
  if (is.null(var_names)) var_names <- paste0("X", seq_len(ncol(X_true)))
  
  nrmse <- vapply(vars, function(j) {
    e <- diff[miss_mask[, j], j]             # alleen de gemiste cellen van kolom j
    sqrt(mean(e^2)) / sd(X_true[, j])
  }, numeric(1))
  
  bias <- vapply(vars, function(j) {
    mean(diff[miss_mask[, j], j])
  }, numeric(1))
  
  data.frame(variable = var_names[vars], nrmse = nrmse, bias = bias,
             row.names = NULL)
}

## Baseline + sanity check: gemiddelde-imputatie ----
##   Verwachting onder MCAR: NRMSE ~ 1, bias ~ 0
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
n         <- 500
p         <- 10
k         <- 5
prop      <- 0.2
strength  <- 2
miss_vars <- 1:3      # could also be: c(1, 2, 7) om ook een ruisvariabele te imputeren in B
driver    <- 4        # buiten miss_vars; in B informatief, in A maakt het niet uit
R         <- 100        
base_seed <- 2026


## Data per scenario ----
gen_data <- function(scenario) {
  switch(scenario,
         A = generate_data_A(n, p),
         B = generate_data_B(n, p))
}

## Methoden: allemaal dezelfde output -> list(X_hat, weights) ----
methods <- list(
  hybrid     = function(Xm) hybrid_kNN(Xm, k = k, weighted = TRUE),
  unweighted = function(Xm) hybrid_kNN(Xm, k = k, weighted = FALSE),
  mean       = function(Xm) list(X_hat = mean_impute(Xm), weights = NULL)
)

## RF-gewichten (list) omzetten naar een tabel, voor de plot in 2.5 ----
weights_to_df <- function(w) {
  do.call(rbind, lapply(names(w), function(t) {
    data.frame(target = t, predictor = names(w[[t]]),
               weight = unname(w[[t]]))
  }))
}

## Eén replicatie ----
one_run <- function(scenario, mechanism, rep) {
  
  set.seed(base_seed + rep)   # zelfde seed per rep -> zelfde complete data over mechanismen
  
  X      <- gen_data(scenario)
  X_miss <- make_missing(X, miss_vars, prop, mechanism,
                         driver = driver, strength = strength)
  
  scores  <- list()
  weights <- NULL
  
  for (m in names(methods)) {
    out   <- methods[[m]](X_miss)            # alle methoden op DEZELFDE X_miss
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

## Grid: 2 scenario's x 3 mechanismen x R replicaties ----
grid <- expand.grid(scenario  = c("A", "B"),
                    mechanism = c("MCAR", "MAR", "MNAR"),
                    rep       = seq_len(R),
                    stringsAsFactors = FALSE)

## Draaien ----
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


## Sanity checks (zet run_checks op TRUE om te draaien) ----
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
  
  # evaluate: mean-imputatie onder MCAR moet NRMSE ~ 1 en bias ~ 0 geven
  print(evaluate(X_B, mean_impute(X_mcar), X_mcar))
  
  # weighted = FALSE moet gelijke gewichten geven
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








