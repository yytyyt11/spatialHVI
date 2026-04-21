# ======================= batch_generate_knn_3groups_halves.R ======================= #
# 三组目录（每组若干 SNP*.csv）：例如 output_SNP50 / output_SNP200 / output_SNP1000
# 为每个 SNPi.csv 生成：
#   simulated_phenotypes.csv, random_effects.csv, residuals.csv, A_used.csv
#
# 统一使用 HALVES 顺序：
#   (L1,...,LJ, R1,...,RJ)
#
# A = kNN(无自环) + 谱归一化；同组所有 replicate 复用同一张 A
# sigma_g^2 以“ROI 对”为单位指定：pair j 的左右半球共用一项
# --------------------------------------------------------------------------- #

suppressPackageStartupMessages({
  library(mvtnorm)
  library(tools)   # file_path_sans_ext
})

# ---------------- Helpers: kNN graph + normalization ---------------- #
make_A_knn <- function(J, k = 2, sigma = 1.0, seed = 123) {
  set.seed(seed)
  if (J <= 0) stop("J must be >= 1")
  if (J == 1) {
    return(matrix(0, 1, 1))
  }
  
  Z  <- matrix(rnorm(J * 2), J, 2)   # 随机 2D 坐标
  D2 <- as.matrix(dist(Z))^2
  W  <- exp(-D2 / (2 * sigma^2))
  diag(W) <- 0
  W  <- 0.5 * (W + t(W))
  
  k_eff <- min(k, max(1, J - 1))
  A <- matrix(0, J, J)
  for (i in 1:J) {
    nn <- order(W[i, ], decreasing = TRUE)[1:k_eff]
    A[i, nn] <- W[i, nn]
  }
  
  A <- pmax(A, t(A))
  diag(A) <- 0
  
  # 连接孤立点
  deg <- rowSums(A)
  iso <- which(deg == 0)
  if (length(iso) > 0) {
    for (i in iso) {
      j <- which.max(W[i, ])
      if (i != j && W[i, j] > 0) A[i, j] <- A[j, i] <- W[i, j]
    }
  }
  
  A
}

normalize_A_spectral <- function(A, eps = 1e-12) {
  A <- 0.5 * (A + t(A))
  diag(A) <- 0
  deg <- rowSums(A)
  Dhi <- diag(1 / sqrt(pmax(deg, eps)), nrow(A))
  B   <- Dhi %*% A %*% Dhi
  lam <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
  lam_max <- max(lam)
  if (!is.finite(lam_max) || lam_max <= 0) lam_max <- 1
  A / lam_max
}

ensure_spd_rho <- function(D_A, A, rho_init = 0.80) {
  spd_ok <- function(r) {
    ev <- eigen(D_A - r * A, symmetric = TRUE, only.values = TRUE)$values
    is.finite(min(ev)) && min(ev) > 0
  }
  
  rho <- min(rho_init, 0.95)
  if (spd_ok(rho)) return(rho)
  
  for (fac in c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)) {
    if (spd_ok(rho * fac)) return(rho * fac)
  }
  
  stop("D_A - rho*A 无法被调整为 SPD，请检查 A。")
}

# SNP 读取（优先空格分隔，退化到逗号）
read_SNP_matrix <- function(path) {
  df <- tryCatch(
    read.csv(path, sep = " ", header = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || ncol(df) == 1) {
    df <- read.csv(path, header = FALSE, check.names = FALSE)
  }
  as.matrix(df)
}

# 构建 K
build_K <- function(SNP_mat) {
  X <- scale(SNP_mat)
  X[!is.finite(X)] <- 0
  L <- ncol(SNP_mat)
  K <- (X %*% t(X)) / L
  K
}

# 将“对级别”的方差向量（长度 J）展开到长度 C=2J，适配 HALVES 顺序：
# (L1,...,LJ,R1,...,RJ)，因此 pair j 对应 trait j 和 j+J
expand_sigma2_pairs_halves <- function(s2_pairs, C) {
  J <- C / 2
  if (length(s2_pairs) != J) {
    warning(sprintf(
      "给定 sigma2_pairs 长度=%d 与 J=%d 不符；将按循环截断/补足。",
      length(s2_pairs), J
    ))
    s2_pairs <- rep(s2_pairs, length.out = J)
  }
  as.numeric(c(s2_pairs, s2_pairs))
}

# 生成一次（针对一个 SNP 文件 + 给定 A/D_A/rho/Sigma/sigma_g^2）
# 统一使用 HALVES 顺序：(L1,...,LJ,R1,...,RJ)
generate_one_halves <- function(SNP_path, out_dir,
                                C, A, D_A, rho, Sigma, sigma2_gc,
                                seed_base = 2025, run_id = 1) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1) Genotype -> K
  G <- read_SNP_matrix(SNP_path)
  n <- nrow(G)
  K <- build_K(G)
  
  # 2) B：每个性状一条 b_c
  set.seed(seed_base + run_id * 11 + C)
  B <- matrix(0, C, n)
  for (c in 1:C) {
    B[c, ] <- as.numeric(rmvnorm(1, sigma = sigma2_gc[c] * K))
  }
  
  # 3) epsilon ~ N(0, Sigma ⊗ (D_A - rho*A)^(-1))  [HALVES]
  spatial_cov <- solve(D_A - rho * A)          # J x J
  full_cov    <- kronecker(Sigma, spatial_cov) # C x C
  
  epsilon <- matrix(0, C, n)
  for (i in 1:n) {
    epsilon[, i] <- as.numeric(rmvnorm(1, sigma = full_cov))
  }
  
  # 4) Y
  Y <- B + epsilon
  
  # 5) 保存
  write.csv(Y,       file.path(out_dir, "simulated_phenotypes.csv"), row.names = FALSE)
  write.csv(B,       file.path(out_dir, "random_effects.csv"),       row.names = FALSE)
  write.csv(epsilon, file.path(out_dir, "residuals.csv"),            row.names = FALSE)
  write.table(A,     file.path(out_dir, "A_used.csv"),
              row.names = FALSE, col.names = FALSE, sep = ",")
  
  invisible(list(n = n, C = C))
}

# ======================= 三组配置（按你的目录改） ======================= #
# 你原来怎么配，这里就怎么配。
# 下面只保留了你上传文件里的一个示例组；如果你有多个组，按同样格式继续加。
group_cfgs <- list(
  output_SNP1000405 = list(
    C = 4,
    rho = 0.00,
    knn_k = 2,
    knn_sigma = 1.0,
    seedA = 1502,
    sigma2_pairs = c(10, 30)
  ),
  output_SNP1000495 = list(
    C = 4,
    rho = 0.90,
    knn_k = 2,
    knn_sigma = 1.0,
    seedA = 1502,
    sigma2_pairs = c(10, 30)
  )
)

# 2x2 残差块协方差（固定，可改）
Sigma_block <- matrix(
  c(10, 5,
    5, 50),
  nrow = 2, byrow = TRUE
)

# ======================= 主循环 ======================= #
summary_rows <- list()
run_counter  <- 0

for (gname in names(group_cfgs)) {
  cfg <- group_cfgs[[gname]]
  C   <- cfg$C
  stopifnot(C %% 2 == 0)
  J   <- C / 2
  
  if (!dir.exists(gname)) {
    warning(sprintf("组目录 '%s' 不存在，跳过。", gname))
    next
  }
  
  # 组内构图（同组复用）
  A_raw <- make_A_knn(J, k = cfg$knn_k, sigma = cfg$knn_sigma, seed = cfg$seedA)
  if (J == 1 && all(A_raw == 0)) A_raw[1, 1] <- 1e-6
  A   <- normalize_A_spectral(A_raw)
  D_A <- diag(rowSums(A))
  rho <- ensure_spd_rho(D_A, A, rho_init = cfg$rho)
  
  # 组内 SNP 列表
  snp_files <- list.files(gname, pattern = "^SNP\\d+\\.csv$", full.names = TRUE)
  if (length(snp_files) == 0) {
    warning(sprintf("组目录 '%s' 下未发现 SNP*.csv，跳过。", gname))
    next
  }
  ord <- order(as.numeric(gsub("^SNP(\\d+)\\.csv$", "\\1", basename(snp_files))))
  snp_files <- snp_files[ord]
  
  # sigma_g^2：从“对级别”展开到 HALVES 顺序
  s2_pairs <- cfg$sigma2_pairs
  sigma2_gc <- expand_sigma2_pairs_halves(s2_pairs, C)
  
  cat(sprintf(
    "\n== 处理组 %s | C=%d (J=%d) | 发现 %d 份 SNP ==\n",
    gname, C, J, length(snp_files)
  ))
  cat(sprintf(
    "A: k=%d, sigma=%.2f, seed=%d | rho(可用)=%.4f | sigma_g^2 pairs = [%s]\n",
    cfg$knn_k, cfg$knn_sigma, cfg$seedA, rho,
    paste(round(s2_pairs, 4), collapse = ", ")
  ))
  
  # 逐 SNP 生成
  for (i in seq_along(snp_files)) {
    snp_path <- snp_files[i]
    base     <- file_path_sans_ext(basename(snp_path))
    out_dir  <- file.path(gname, base)   # e.g. output_SNP50/SNP1/
    
    t0 <- proc.time()[3]
    generate_one_halves(
      SNP_path = snp_path,
      out_dir = out_dir,
      C = C,
      A = A,
      D_A = D_A,
      rho = rho,
      Sigma = Sigma_block,
      sigma2_gc = sigma2_gc,
      seed_base = 777 + which(names(group_cfgs) == gname) * 1000,
      run_id = i
    )
    t1 <- proc.time()[3]
    
    run_counter <- run_counter + 1
    summary_rows[[run_counter]] <- data.frame(
      group = gname,
      snp_file = basename(snp_path),
      out_dir  = out_dir,
      C = C,
      J = J,
      rho_used = rho,
      time_sec = round(t1 - t0, 4),
      stringsAsFactors = FALSE
    )
    
    cat(sprintf("  - 完成 %-12s -> %s (%.3fs)\n",
                basename(snp_path), out_dir, t1 - t0))
  }
}

# 保存总汇
if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  write.csv(summary_df, "batch_generation_summary2.csv", row.names = FALSE)
  cat(sprintf("\n全部完成：共 %d 次生成。概要已写入 batch_generation_summary2.csv\n", nrow(summary_df)))
} else {
  cat("\n没有任何生成（可能未找到组目录或 SNP 文件）。\n")
}
