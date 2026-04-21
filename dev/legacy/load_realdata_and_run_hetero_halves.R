## ====================== load_realdata_and_run_hetero_halves.R ====================== ##
## Usage:
##   source("vi_halves_unified.R")
##   source("load_realdata_and_run_hetero_halves.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

## 1) Read data
BrainMeasures <- as.matrix(read.csv("BrainMeasure_data.csv", header = TRUE, row.names = 1))
load("FreeSurfer_Data.RData")   # should provide FreeSurfer_list_Data

## 2) Genotype matrix: subjects x SNPs
## Here we keep the same slice as your current script
X <- as.matrix(FreeSurfer_list_Data[1:632, 2:487])
M <- nrow(X)   # number of subjects = 632
L <- ncol(X)   # number of SNPs = 486

## 3) Phenotypes: row-wise z-score, then use C x M matrix
Y_raw <- BrainMeasures[, 1:M, drop = FALSE]
Y <- t(apply(Y_raw, 1, scale))
Y[!is.finite(Y)] <- 0

C <- nrow(Y)
stopifnot(C %% 2 == 0)
J <- C / 2

## 4) Build K
std_X <- safe_scale_cols(X)
K <- (std_X %*% t(std_X)) / L

## 5) Build A under halves ordering using the same rule as the core file
A <- build_A_from_Y_halves(Y, k = 2)
D_A <- diag(rowSums(A))  # for inspection if needed

## 6) Run heteroscedastic VI (real-data entry)
res <- run_vi_hetero_from_mats(
  Y = Y,
  K = K,
  A_mat = A,
  rho_grid = seq(0.55, 0.99, by = 0.01),
  max_iter = 10000,
  tol = 1e-4,
  verbose = TRUE
)

## 7) ELBO trace
plot(
  seq_along(res$elbo), res$elbo, type = "l",
  main = "ELBO (Heteroscedastic, halves)",
  xlab = "Iteration", ylab = "ELBO"
)

## 8) Heritability
h <- heritability_from_vi(res)

cat(sprintf(
  "\n[Summary]\n  E[rho] = %.4f\n  Broad H^2 = %.4f\n  mean(narrow h^2) = %.4f (sd=%.4f)\n",
  res$E_rho, h$H2, mean(h$h2), sd(h$h2)
))
cat("  First 10 narrow h^2:\n")
print(round(head(h$h2, 10), 4))

## 9) Save posterior summaries
save_vi_result(res, prefix = "brain_hetero_halves")