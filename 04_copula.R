# =============================================================
# 04_copula.R
# Bivariate and trivariate copula fitting, tail dependence
# coefficients, R-vine structure selection, and PIT scatter
# plots.
#
# Requires: pit_uniforms.csv produced by 03_garch.R
#           OR the U matrix still in your R session
# Output:   bivariate_copula_results.csv,
#           trivariate_comparison.csv,
#           pit_scatter_pairs.png,
#           copula_results.RData
# =============================================================

library(copula)
library(VineCopula)

# --- Load PIT uniforms if starting a fresh session ---
# U <- as.matrix(read.csv("pit_uniforms.csv"))

# =============================================================
# PART 1: Bivariate copula fitting
# =============================================================

pairs_list <- list(
  c("BTC",   "ETH"),
  c("BTC",   "SP500"),
  c("ETH",   "SP500")
)

biv_results <- data.frame()

for (pair in pairs_list) {
  u_pair    <- U[, pair]
  pair_name <- paste(pair[1], pair[2], sep = "-")
  n         <- nrow(u_pair)

  gauss_fit <- fitCopula(normalCopula(dim = 2),  u_pair, method = "ml")
  t_fit     <- fitCopula(tCopula(dim = 2),       u_pair, method = "ml")
  clay_fit  <- fitCopula(claytonCopula(dim = 2), u_pair, method = "ml")
  gumb_fit  <- fitCopula(gumbelCopula(dim = 2),  u_pair, method = "ml")

  aic_fn <- function(f) -2 * f@loglik + 2 * length(coef(f))
  bic_fn <- function(f) -2 * f@loglik + log(n) * length(coef(f))

  biv_results <- rbind(biv_results, data.frame(
    Pair   = pair_name,
    Copula = c("Gaussian", "Student-t", "Clayton", "Gumbel"),
    Param1 = c(coef(gauss_fit)[1], coef(t_fit)[1],
               coef(clay_fit)[1], coef(gumb_fit)[1]),
    Param2 = c(NA, coef(t_fit)[2], NA, NA),
    LogLik = c(gauss_fit@loglik, t_fit@loglik,
               clay_fit@loglik, gumb_fit@loglik),
    AIC    = c(aic_fn(gauss_fit), aic_fn(t_fit),
               aic_fn(clay_fit), aic_fn(gumb_fit)),
    BIC    = c(bic_fn(gauss_fit), bic_fn(t_fit),
               bic_fn(clay_fit), bic_fn(gumb_fit))
  ))
}

print(biv_results, digits = 4, row.names = FALSE)
write.csv(biv_results, "bivariate_copula_results.csv",
          row.names = FALSE)

# =============================================================
# PART 2: Tail dependence coefficients
# =============================================================

cat("\n===== TAIL DEPENDENCE COEFFICIENTS =====\n")
for (pair in pairs_list) {
  u_pair    <- U[, pair]
  pair_name <- paste(pair[1], pair[2], sep = "-")

  t2  <- fitCopula(tCopula(dim = 2),       u_pair, method = "ml")
  cl2 <- fitCopula(claytonCopula(dim = 2), u_pair, method = "ml")
  gu2 <- fitCopula(gumbelCopula(dim = 2),  u_pair, method = "ml")

  rho    <- coef(t2)[1]
  df     <- coef(t2)[2]
  lam_t  <- 2 * pt(-sqrt((df + 1) * (1 - rho) / (1 + rho)), df = df + 1)
  lam_cl <- 2^(-1 / coef(cl2)[1])
  lam_gu <- 2 - 2^(1 / coef(gu2)[1])

  cat(pair_name,
      ": t lambda (L=U) =", round(lam_t,  4),
      " | Clayton lambdaL =", round(lam_cl, 4),
      " | Gumbel lambdaU =",  round(lam_gu, 4), "\n")
}

# =============================================================
# PART 3: Trivariate model comparison
# =============================================================

cat("\n===== TRIVARIATE COPULA COMPARISON =====\n")

n <- nrow(U)

# 3D Gaussian copula
gauss3_fit <- fitCopula(
  normalCopula(dim = 3, dispstr = "un"), U, method = "ml"
)

# 3D Student-t copula
t3_fit <- fitCopula(
  tCopula(dim = 3, dispstr = "un"), U, method = "ml"
)

# R-vine copula (automated structure and family selection)
vine_fit <- RVineStructureSelect(
  U,
  type          = 0,
  familyset     = NA,
  selectioncrit = "AIC"
)
summary(vine_fit)

vine_loglik <- as.numeric(RVineLogLik(U, vine_fit)$loglik)
vine_npar   <- sum(vine_fit$family  > 0) +
               sum(vine_fit$family2 > 0)

comparison <- data.frame(
  Model      = c("3D Gaussian", "3D Student-t", "R-Vine"),
  Parameters = c(length(coef(gauss3_fit)),
                 length(coef(t3_fit)),
                 vine_npar),
  LogLik     = c(gauss3_fit@loglik,
                 t3_fit@loglik,
                 vine_loglik),
  AIC = c(
    -2 * gauss3_fit@loglik + 2 * length(coef(gauss3_fit)),
    -2 * t3_fit@loglik     + 2 * length(coef(t3_fit)),
    -2 * vine_loglik       + 2 * vine_npar
  ),
  BIC = c(
    -2 * gauss3_fit@loglik + log(n) * length(coef(gauss3_fit)),
    -2 * t3_fit@loglik     + log(n) * length(coef(t3_fit)),
    -2 * vine_loglik       + log(n) * vine_npar
  )
)
print(comparison, digits = 4, row.names = FALSE)
write.csv(comparison, "trivariate_comparison.csv",
          row.names = FALSE)

# =============================================================
# PART 4: Figure 7 - Pairwise PIT scatter plots
# =============================================================

png("pit_scatter_pairs.png", width = 1600, height = 1600, res = 200)
pairs(
  U,
  labels = c("BTC", "ETH", "S&P 500"),
  main   = "Pairwise PIT Scatter Plots",
  pch    = ".",
  col    = rgb(0, 0, 0, 0.3)
)
dev.off()

# --- Save copula objects for use in 05_risk.R ---
save(vine_fit, gauss3_fit, t3_fit, U,
     file = "copula_results.RData")

cat("Done. Saved: bivariate_copula_results.csv,",
    "trivariate_comparison.csv, pit_scatter_pairs.png,",
    "copula_results.RData\n")
