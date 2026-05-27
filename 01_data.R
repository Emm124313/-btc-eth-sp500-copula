# =============================================================
# 01_data.R
# Data retrieval and log-return construction
#
# Downloads daily closing prices for BTC, ETH, and S&P 500
# from Yahoo Finance, aligns on common trading dates, and
# computes continuously compounded log returns.
#
# Output files:
#   aligned_prices.csv  -- 1,277 rows of aligned closing prices
#   log_returns.csv     -- 1,276 rows of daily log returns
# =============================================================

library(quantmod)

# --- Download raw price series ---
getSymbols("BTC-USD", src = "yahoo",
           from = "2020-01-01", to = "2025-01-31")
getSymbols("ETH-USD", src = "yahoo",
           from = "2020-01-01", to = "2025-01-31")
getSymbols("^GSPC",   src = "yahoo",
           from = "2020-01-01", to = "2025-01-31")

# --- Extract adjusted closing prices ---
btc <- Ad(`BTC-USD`)
eth <- Ad(`ETH-USD`)
sp  <- Ad(GSPC)

# --- Align on shared calendar dates (inner join) ---
prices <- merge(btc, eth, sp, join = "inner")
colnames(prices) <- c("BTC", "ETH", "SP500")
cat("Aligned price observations:", nrow(prices), "\n")  # expected: 1,277

# --- Compute continuously compounded log returns ---
returns <- diff(log(prices))[-1, ]
colnames(returns) <- c("BTC", "ETH", "SP500")
cat("Return observations:", nrow(returns), "\n")         # expected: 1,276

# --- Save to CSV ---
write.csv(
  data.frame(Date = index(prices), coredata(prices)),
  "aligned_prices.csv",
  row.names = FALSE
)
write.csv(
  data.frame(Date = index(returns), coredata(returns)),
  "log_returns.csv",
  row.names = FALSE
)

cat("Files saved: aligned_prices.csv, log_returns.csv\n")
