# =============================================================
# 03_garch.R
# GARCH(1,1) marginal model estimation with Hansen skewed
# Student's t innovations, PIT diagnostics, and figures.
#
# Requires: 'returns' object from 01_data.R
#           OR load log_returns.csv (see commented lines below)
# Output:   pit_uniforms.csv, conditional_volatility.png,
#           garch_qqplot.png, pit_histograms.png
# =============================================================

library(rugarch)

# --- Load returns if starting a fresh session ---
# library(xts)
# returns_df <- read.csv("log_returns.csv")
# returns    <- xts(returns_df[,-1],
#                   order.by = as.Date(returns_df$Date))

# --- GARCH(1,1) specification with skewed-t innovations ---
garch_spec <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "sstd"
)

# --- Fit to each series ---
fit_btc <- ugarchfit(spec = garch_spec, data = returns[, "BTC"])
fit_eth <- ugarchfit(spec = garch_spec, data = returns[, "ETH"])
fit_sp  <- ugarchfit(spec = garch_spec, data = returns[, "SP500"])

show(fit_btc)
show(fit_eth)
show(fit_sp)

# --- Parameter comparison table ---
extract_params <- function(fit, name) {
  cf <- coef(fit)
  data.frame(
    Asset       = name,
    mu          = cf["mu"],
    omega       = cf["omega"],
    alpha1      = cf["alpha1"],
    beta1       = cf["beta1"],
    skew        = cf["skew"],
    shape       = cf["shape"],
    persistence = cf["alpha1"] + cf["beta1"],
    LogLik      = likelihood(fit),
    AIC         = infocriteria(fit)[1],
    BIC         = infocriteria(fit)[2]
  )
}
params_table <- rbind(
  extract_params(fit_btc, "BTC"),
  extract_params(fit_eth, "ETH"),
  extract_params(fit_sp,  "SP500")
)
rownames(params_table) <- NULL
print(params_table, digits = 6)

# --- Probability Integral Transform (PIT) ---
pit_btc <- pit(fit_btc)
pit_eth <- pit(fit_eth)
pit_sp  <- pit(fit_sp)

U <- cbind(
  BTC   = as.numeric(pit_btc),
  ETH   = as.numeric(pit_eth),
  SP500 = as.numeric(pit_sp)
)
cat("PIT value ranges (should all be in (0,1)):\n")
print(apply(U, 2, range))

# --- KS uniformity tests on PIT values ---
cat("\n===== KS UNIFORMITY TEST =====\n")
for (col in colnames(U)) {
  res <- ks.test(U[, col], "punif", 0, 1)
  cat(col, ": D =", round(res$statistic, 4),
      ", p =", format.pval(res$p.value, digits = 4), "\n")
}

# --- Ljung-Box on PIT values and squared PIT values ---
cat("\n===== LJUNG-BOX ON PIT =====\n")
for (col in colnames(U)) {
  r1 <- Box.test(U[, col],    lag = 20, type = "Ljung-Box")
  r2 <- Box.test(U[, col]^2, lag = 20, type = "Ljung-Box")
  cat(col, "level p =", format.pval(r1$p.value, digits = 4),
      " | squared p =", format.pval(r2$p.value, digits = 4), "\n")
}

# --- Save PIT uniforms for copula stage ---
write.csv(U, "pit_uniforms.csv", row.names = FALSE)

# --- Figure 4: conditional volatility ---
png("conditional_volatility.png", width = 1800, height = 1200, res = 150)
par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))
plot(sigma(fit_btc),
     main = "Conditional Volatility: BTC",
     ylab = expression(sigma[t]))
plot(sigma(fit_eth),
     main = "Conditional Volatility: ETH",
     ylab = expression(sigma[t]))
plot(sigma(fit_sp),
     main = "Conditional Volatility: S&P 500",
     ylab = expression(sigma[t]))
dev.off()

# --- Figure 5: QQ plots of standardised GARCH residuals ---
png("garch_qqplot.png", width = 1800, height = 600, res = 150)
par(mfrow = c(1, 3))
plot(fit_btc, which = 9)
plot(fit_eth, which = 9)
plot(fit_sp,  which = 9)
dev.off()

# --- Figure 6: PIT diagnostic histograms ---
png("pit_histograms.png", width = 1800, height = 600, res = 150)
par(mfrow = c(1, 3))
hist(U[, "BTC"],   breaks = 30, main = "PIT: BTC",
     xlab = "u", col = "lightblue",   freq = FALSE)
abline(h = 1, col = "red", lwd = 2)
hist(U[, "ETH"],   breaks = 30, main = "PIT: ETH",
     xlab = "u", col = "lightgreen",  freq = FALSE)
abline(h = 1, col = "red", lwd = 2)
hist(U[, "SP500"], breaks = 30, main = "PIT: S&P 500",
     xlab = "u", col = "lightyellow", freq = FALSE)
abline(h = 1, col = "red", lwd = 2)
dev.off()

# --- Save fitted GARCH objects for use in 05_risk.R ---
save(fit_btc, fit_eth, fit_sp, garch_spec, U,
     file = "garch_fits.RData")

cat("Done. Saved: pit_uniforms.csv, garch_fits.RData,",
    "conditional_volatility.png, garch_qqplot.png,",
    "pit_histograms.png\n")
