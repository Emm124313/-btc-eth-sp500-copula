# =============================================================
# 02_diagnostics.R
# Preliminary analysis: descriptive statistics, normality,
# stationarity, serial dependence, ARCH-LM tests, and
# all preliminary figures.
#
# Requires: log_returns.csv produced by 01_data.R
#           OR the 'returns' object still in your R session.
# =============================================================

library(moments)
library(tseries)
library(FinTS)

# --- Load returns if starting a fresh session ---
# returns_df <- read.csv("log_returns.csv")
# returns    <- xts(returns_df[,-1],
#                   order.by = as.Date(returns_df$Date))

# --- Descriptive statistics ---
stats_table <- data.frame(
  Asset       = c("BTC", "ETH", "S&P 500"),
  Obs         = nrow(returns),
  Mean        = round(colMeans(returns),        6),
  Std_Dev     = round(apply(returns, 2, sd),    6),
  Min         = round(apply(returns, 2, min),   6),
  Max         = round(apply(returns, 2, max),   6),
  Skewness    = round(apply(returns, 2, skewness),          4),
  Ex_Kurtosis = round(apply(returns, 2, kurtosis) - 3,      4)
)
print(stats_table)
write.csv(stats_table, "descriptive_statistics.csv",
          row.names = FALSE)

# --- Jarque-Bera normality tests ---
cat("\n===== JARQUE-BERA TEST =====\n")
print(jarque.test(as.numeric(returns$BTC)))
print(jarque.test(as.numeric(returns$ETH)))
print(jarque.test(as.numeric(returns$SP500)))

# --- ADF stationarity tests ---
cat("\n===== ADF STATIONARITY TEST =====\n")
print(adf.test(as.numeric(returns$BTC)))
print(adf.test(as.numeric(returns$ETH)))
print(adf.test(as.numeric(returns$SP500)))

# --- Ljung-Box on returns ---
cat("\n===== LJUNG-BOX TEST (returns) =====\n")
print(Box.test(as.numeric(returns$BTC),   lag = 20, type = "Ljung-Box"))
print(Box.test(as.numeric(returns$ETH),   lag = 20, type = "Ljung-Box"))
print(Box.test(as.numeric(returns$SP500), lag = 20, type = "Ljung-Box"))

# --- Ljung-Box on squared returns (volatility clustering) ---
cat("\n===== LJUNG-BOX TEST (squared returns) =====\n")
print(Box.test(as.numeric(returns$BTC)^2,   lag = 20, type = "Ljung-Box"))
print(Box.test(as.numeric(returns$ETH)^2,   lag = 20, type = "Ljung-Box"))
print(Box.test(as.numeric(returns$SP500)^2, lag = 20, type = "Ljung-Box"))

# --- ARCH-LM test ---
cat("\n===== ARCH-LM TEST (12 lags) =====\n")
print(ArchTest(as.numeric(returns$BTC),   lags = 12))
print(ArchTest(as.numeric(returns$ETH),   lags = 12))
print(ArchTest(as.numeric(returns$SP500), lags = 12))

# --- Correlation matrices ---
cat("\n===== PEARSON CORRELATION =====\n")
print(round(cor(returns, method = "pearson"),  4))
cat("\n===== KENDALL TAU =====\n")
print(round(cor(returns, method = "kendall"),  4))
cat("\n===== SPEARMAN RHO =====\n")
print(round(cor(returns, method = "spearman"), 4))

# --- Figure 1: time-series of log returns ---
png("returns_timeseries.png", width = 1800, height = 1200, res = 150)
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
plot(index(returns), as.numeric(returns$BTC),
     type = "l", col = "blue",
     main = "BTC Daily Log Returns",
     xlab = "Date", ylab = "Log Return")
abline(h = 0, col = "red", lty = 2)
plot(index(returns), as.numeric(returns$ETH),
     type = "l", col = "darkgreen",
     main = "ETH Daily Log Returns",
     xlab = "Date", ylab = "Log Return")
abline(h = 0, col = "red", lty = 2)
plot(index(returns), as.numeric(returns$SP500),
     type = "l", col = "black",
     main = "S&P 500 Daily Log Returns",
     xlab = "Date", ylab = "Log Return")
abline(h = 0, col = "red", lty = 2)
dev.off()

# --- Figure 2: histograms with normal density overlay ---
png("return_histograms.png", width = 1800, height = 600, res = 150)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
hist(as.numeric(returns$BTC), breaks = 50, probability = TRUE,
     main = "BTC Returns", xlab = "Log Return", col = "lightblue")
curve(dnorm(x, mean = mean(returns$BTC), sd = sd(returns$BTC)),
      add = TRUE, col = "red", lwd = 2)
hist(as.numeric(returns$ETH), breaks = 50, probability = TRUE,
     main = "ETH Returns", xlab = "Log Return", col = "lightgreen")
curve(dnorm(x, mean = mean(returns$ETH), sd = sd(returns$ETH)),
      add = TRUE, col = "red", lwd = 2)
hist(as.numeric(returns$SP500), breaks = 50, probability = TRUE,
     main = "S&P 500 Returns", xlab = "Log Return", col = "lightgray")
curve(dnorm(x, mean = mean(returns$SP500), sd = sd(returns$SP500)),
      add = TRUE, col = "red", lwd = 2)
dev.off()

# --- Figure 3: normal QQ plots ---
png("qqplot_returns.png", width = 1800, height = 600, res = 150)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
qqnorm(as.numeric(returns$BTC),   main = "QQ Plot: BTC")
qqline(as.numeric(returns$BTC),   col = "red", lwd = 2)
qqnorm(as.numeric(returns$ETH),   main = "QQ Plot: ETH")
qqline(as.numeric(returns$ETH),   col = "red", lwd = 2)
qqnorm(as.numeric(returns$SP500), main = "QQ Plot: S&P 500")
qqline(as.numeric(returns$SP500), col = "red", lwd = 2)
dev.off()

cat("Figures saved: returns_timeseries.png,",
    "return_histograms.png, qqplot_returns.png\n")
