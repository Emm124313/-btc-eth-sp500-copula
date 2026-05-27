# =============================================================
# 05_risk.R
# Rolling Monte Carlo VaR and ES estimation, followed by
# Kupiec unconditional coverage and Christoffersen conditional
# coverage backtests.
#
# Requires: garch_fits.RData from 03_garch.R
#           log_returns.csv  from 01_data.R
# Output:   var_es_results.csv, backtest_results.csv,
#           rolling_var_99.png
# =============================================================

library(rugarch)
library(copula)
library(VineCopula)

# --- Load fitted GARCH objects and return data ---
load("garch_fits.RData")   # loads: fit_btc, fit_eth, fit_sp,
                            #        garch_spec, U
returns_df <- read.csv("log_returns.csv")
R_mat <- as.matrix(returns_df[, c("BTC", "ETH", "SP500")])

# =============================================================
# PART 1: Rolling Monte Carlo simulation
# =============================================================

set.seed(42)
N_sim       <- 10000
T1          <- 250
alpha_vec   <- c(0.95, 0.99)
T_obs       <- nrow(R_mat)
model_names <- c("Gaussian", "t", "RVine")
nmodels     <- length(model_names)

VaR_all <- array(
  NA,
  dim      = c(T1, nmodels, length(alpha_vec)),
  dimnames = list(NULL, model_names,
                  paste0(alpha_vec * 100, "%"))
)
ES_all     <- VaR_all
violations <- matrix(
  0L,
  nrow     = nmodels,
  ncol     = length(alpha_vec),
  dimnames = list(model_names,
                  paste0(alpha_vec * 100, "%"))
)

for (day in 1:T1) {

  t_end <- T_obs - T1 + day - 1
  t_new <- t_end + 1
  R_in  <- R_mat[1:t_end, ]

  # Re-fit GARCH marginals on expanding in-sample window
  f1 <- ugarchfit(garch_spec, R_in[, 1], solver = "hybrid")
  f2 <- ugarchfit(garch_spec, R_in[, 2], solver = "hybrid")
  f3 <- ugarchfit(garch_spec, R_in[, 3], solver = "hybrid")

  # One-step-ahead conditional volatility forecasts
  sigma_hat <- c(
    as.numeric(sigma(ugarchforecast(f1, n.ahead = 1))),
    as.numeric(sigma(ugarchforecast(f2, n.ahead = 1))),
    as.numeric(sigma(ugarchforecast(f3, n.ahead = 1)))
  )
  mu_hat <- c(coef(f1)["mu"], coef(f2)["mu"], coef(f3)["mu"])

  # PIT helper functions
  pit_fn <- function(z, f) {
    cf <- coef(f)
    as.numeric(pdist("sstd", q = z, mu = 0, sigma = 1,
                     skew = cf["skew"], shape = cf["shape"]))
  }
  qpit_fn <- function(p, f) {
    cf <- coef(f)
    as.numeric(qdist("sstd", p = p, mu = 0, sigma = 1,
                     skew = cf["skew"], shape = cf["shape"]))
  }

  # Build in-sample PIT uniforms
  U_in <- cbind(
    pit_fn(as.numeric(residuals(f1, standardize = TRUE)), f1),
    pit_fn(as.numeric(residuals(f2, standardize = TRUE)), f2),
    pit_fn(as.numeric(residuals(f3, standardize = TRUE)), f3)
  )
  U_in <- pmax(pmin(U_in, 1 - 1e-6), 1e-6)  # numerical safety

  # Fit three copula models
  cop_g  <- tryCatch(
    fitCopula(normalCopula(dim = 3, dispstr = "un"),
              U_in, method = "ml"),
    error = function(e) NULL
  )
  cop_t  <- tryCatch(
    fitCopula(tCopula(dim = 3, dispstr = "un"),
              U_in, method = "ml"),
    error = function(e) NULL
  )
  cop_rv <- tryCatch(
    RVineStructureSelect(U_in, familyset = NA,
                         selectioncrit = "AIC"),
    error = function(e) NULL
  )

  models <- list(cop_g, cop_t, cop_rv)

  for (m in seq_along(models)) {
    mod <- models[[m]]
    if (is.null(mod)) next

    # Simulate N_sim joint uniform vectors
    U_sim <- if (m <= 2) {
      rCopula(N_sim, mod@copula)
    } else {
      RVineSim(N_sim, mod)
    }

    # Invert PIT to get simulated standardised innovations
    z_sim <- cbind(
      qpit_fn(U_sim[, 1], f1),
      qpit_fn(U_sim[, 2], f2),
      qpit_fn(U_sim[, 3], f3)
    )

    # Construct simulated one-day returns
    r_sim <- sweep(
      sweep(z_sim, 2, sigma_hat, "*"),
      2, mu_hat, "+"
    )

    # Equally weighted portfolio loss
    loss <- -rowMeans(r_sim)

    # VaR and ES at each confidence level
    for (a in seq_along(alpha_vec)) {
      VaR_all[day, m, a] <- quantile(loss, probs = alpha_vec[a])
      ES_all[day, m, a]  <- mean(
        loss[loss > VaR_all[day, m, a]]
      )
    }
  }

  # Realised portfolio loss
  loss_actual <- -mean(R_mat[t_new, ])

  # Count violations
  for (m in 1:nmodels)
    for (a in seq_along(alpha_vec))
      if (!is.na(VaR_all[day, m, a]) &&
          loss_actual > VaR_all[day, m, a])
        violations[m, a] <- violations[m, a] + 1L

  if (day %% 50 == 0)
    cat(sprintf("Day %d / %d complete\n", day, T1))
}

# =============================================================
# PART 2: Average VaR and ES table
# =============================================================

cat("\n===== AVERAGE VaR AND ES =====\n")
var_es_rows <- data.frame()
for (m in 1:nmodels) {
  for (a in seq_along(alpha_vec)) {
    var_es_rows <- rbind(var_es_rows, data.frame(
      Model      = model_names[m],
      Level      = paste0(alpha_vec[a] * 100, "%"),
      Avg_VaR    = round(mean(VaR_all[, m, a], na.rm = TRUE), 5),
      Avg_ES     = round(mean(ES_all[,  m, a], na.rm = TRUE), 5),
      Violations = violations[m, a]
    ))
  }
}
print(var_es_rows, row.names = FALSE)
write.csv(var_es_rows, "var_es_results.csv", row.names = FALSE)

# =============================================================
# PART 3: Kupiec unconditional coverage test
# =============================================================

kupiec_test <- function(x, T1, p) {
  x_adj <- max(x, 1e-4)
  lr    <- -2 * (
    x       * log(p)       + (T1 - x) * log(1 - p) -
    x_adj   * log(x_adj / T1) -
    (T1 - x) * log(1 - x_adj / T1)
  )
  pval <- 1 - pchisq(lr, df = 1)
  list(LR = round(lr, 4), pval = round(pval, 4),
       result = ifelse(pval < 0.05, "REJECT", "Pass"))
}

cat("\n===== KUPIEC UNCONDITIONAL COVERAGE TEST =====\n")
backtest_rows <- data.frame()
for (m in 1:nmodels) {
  for (a in seq_along(alpha_vec)) {
    res <- kupiec_test(violations[m, a], T1, 1 - alpha_vec[a])
    cat(sprintf(
      "%-10s %s%%  violations=%d  rate=%.4f  LR=%.4f  p=%.4f  %s\n",
      model_names[m], alpha_vec[a] * 100,
      violations[m, a], violations[m, a] / T1,
      res$LR, res$pval, res$result
    ))
    backtest_rows <- rbind(backtest_rows, data.frame(
      Model      = model_names[m],
      Level      = paste0(alpha_vec[a] * 100, "%"),
      Expected   = T1 * (1 - alpha_vec[a]),
      Observed   = violations[m, a],
      Rate       = round(violations[m, a] / T1, 4),
      LR_uc      = res$LR,
      p_value    = res$pval,
      Result     = res$result
    ))
  }
}
write.csv(backtest_rows, "backtest_results.csv", row.names = FALSE)

# =============================================================
# PART 4: Figure 8 - Rolling 99% VaR vs realised losses
# =============================================================

returns_vec <- -rowMeans(
  R_mat[(T_obs - T1 + 1):T_obs, ]
)

png("rolling_var_99.png", width = 1800, height = 800, res = 150)
par(mar = c(4, 4, 3, 1))
plot(
  1:T1, returns_vec,
  type = "l", col = "black", lwd = 1,
  main = "Rolling 99% VaR Forecasts vs Realised Portfolio Losses",
  xlab = "Out-of-sample day",
  ylab = "Portfolio loss"
)
lines(1:T1, VaR_all[, "Gaussian", "99%"],
      col = "gray50", lwd = 2)
lines(1:T1, VaR_all[, "t",        "99%"],
      col = "blue",   lwd = 2)
lines(1:T1, VaR_all[, "RVine",    "99%"],
      col = "red",    lwd = 2)
legend(
  "topright",
  legend = c("Realised loss", "Gaussian",
             "Student-t", "R-Vine"),
  col    = c("black", "gray50", "blue", "red"),
  lwd    = 2, bty = "n"
)
dev.off()

cat("Done. Saved: var_es_results.csv, backtest_results.csv,",
    "rolling_var_99.png\n")

