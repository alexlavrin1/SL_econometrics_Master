# ==============================================================================
# Part 2 — Assignment Statistical Learning
# Author: Sebastiaan
# ==============================================================================

# Libraries ---------------------------------------------------------------

library(MASS)


# 2.1 Data generating processes -------------------------------------------

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

## Scenario B: signaalblok + ruis ----
##   X1..X_signal = Z + e_i  (gedeelde latente factor)
##   X_(signal+1)..Xp        = onafhankelijke ruis
##   Impliciete correlatie binnen het blok: Var(Z) / (Var(Z) + var_e)
##   Verwachting: RF legt bijna al het gewicht op het signaalblok,
##   ongewogen afstand wordt verstoord door de ruisdimensies.

generate_data_B <- function(n, p = 10, n_signal = 4, var_e = 0.3) {
  Z <- rnorm(n, mean = 0, sd = 1)
  
  signal <- sapply(seq_len(n_signal), function(i) {
    Z + rnorm(n, mean = 0, sd = sqrt(var_e))
  })
  
  noise <- matrix(rnorm(n * (p - n_signal), mean = 0, sd = 1),
                  nrow = n, ncol = p - n_signal)
  
  X <- cbind(signal, noise)
  colnames(X) <- paste0("X", seq_len(p))
  X
}


# 2.2 Missing data mechanisms ---------------------------------------------

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


# Test --------------------------------------------------------------------

set.seed(2026)

X_A <- generate_data_A(n = 500, p = 10)
X_B <- generate_data_B(n = 500, p = 10)

round(cor(X_B), 2)   # blok X1-X4 hoog gecorreleerd, rest ~ 0

miss_vars <- 1:3
driver    <- 4       # informatieve variabele, buiten miss_vars

X_mcar <- make_missing(X_B, miss_vars, 0.2, "MCAR")
X_mar  <- make_missing(X_B, miss_vars, 0.2, "MAR",  driver = driver, strength = 2)
X_mnar <- make_missing(X_B, miss_vars, 0.2, "MNAR", driver = driver, strength = 2)

cat("\n--- MCAR ---\n"); check_missing(X_B, X_mcar, miss_vars, driver)
cat("\n--- MAR  ---\n"); check_missing(X_B, X_mar,  miss_vars, driver)
cat("\n--- MNAR ---\n"); check_missing(X_B, X_mnar, miss_vars, driver)



# 2.3 Evaluation metric ---------------------------------------------------
# Wacht op part 1: aggregatiemethode bepaalt metric (mean -> RMSE, median -> MAE)





# 2.4 Simulation kNN ------------------------------------------------------






# 2.5 Summarize data  ------------------------------------------------------






