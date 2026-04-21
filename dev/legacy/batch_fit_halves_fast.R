source("vi_halves_unified_fast.R")

group_dir <- "output_SNP1000495"
resume_from_snp <- 1L
skip_completed <- FALSE

snp_files <- list.files(group_dir, pattern = "^SNP\\d+\\.csv$", full.names = TRUE)
ord <- order(as.numeric(gsub("^SNP(\\d+)\\.csv$", "\\1", basename(snp_files))))
snp_files <- snp_files[ord]
snp_ids <- as.integer(gsub("^SNP(\\d+)\\.csv$", "\\1", basename(snp_files)))

build_summary_row <- function(snp_name, res, a_path) {
  A_true <- as.matrix(read.csv(a_path, header = FALSE))
  h <- heritability_from_vi(res, A_mat = A_true)

  data.frame(
    snp = snp_name,
    E_rho = res$E_rho,
    H2 = h$H2,
    mean_h2 = mean(h$h2),
    sd_h2 = sd(h$h2)
  )
}

summary_list <- list()

for (i in seq_along(snp_files)) {
  snp_path <- snp_files[i]
  snp_id <- snp_ids[i]
  snp_name <- tools::file_path_sans_ext(basename(snp_path))
  data_dir <- file.path(group_dir, snp_name)

  y_path <- file.path(data_dir, "simulated_phenotypes.csv")
  a_path <- file.path(data_dir, "A_used.csv")
  result_path <- file.path(data_dir, "vi_halves_fast_result.rds")

  if (!file.exists(y_path) || !file.exists(a_path)) {
    cat("Skipping", snp_name, ": missing input files.\n")
    next
  }

  if (skip_completed && file.exists(result_path)) {
    cat("\n===== Reusing", snp_name, "=====\n")
    res_existing <- readRDS(result_path)
    summary_list[[length(summary_list) + 1]] <- build_summary_row(snp_name, res_existing, a_path)
    next
  }

  if (snp_id < resume_from_snp) {
    cat("\n===== Skipping", snp_name, "(before resume point) =====\n")
    next
  }

  cat("\n===== Running", snp_name, "=====\n")

  Y <- read_numeric_csv_matrix(y_path, row_names = FALSE)

  read_SNP_matrix <- function(path) {
    G <- tryCatch(
      read.table(path, header = FALSE, sep = "", check.names = FALSE),
      error = function(e) NULL
    )
    if (is.null(G) || ncol(G) == 1L) {
      G <- read.csv(path, header = FALSE, check.names = FALSE)
    }
    as.matrix(G)
  }

  G <- read_SNP_matrix(snp_path)
  X <- safe_scale_cols(G)
  L <- ncol(X)
  K <- (X %*% t(X)) / L

  A_true <- as.matrix(read.csv(a_path, header = FALSE))

  fit <- run_sim_from_mats_hetero(
    Y = Y,
    K = K,
    A_true = A_true,
    rho_grid = seq(0.00, 0.99, by = 0.01),
    max_iter = 5000,
    tol = 1e-4,
    verbose = TRUE,
    prefix = file.path(data_dir, "vi_halves_fast"),
    init_h2 = 0.15,
    init_ig_shape = 3.0,
    init_iw_df = 8
  )

  summary_list[[length(summary_list) + 1]] <- build_summary_row(snp_name, fit$res, a_path)
}

if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  summary_df <- summary_df[order(as.integer(gsub("^SNP", "", summary_df$snp))), ]
  write.csv(summary_df, file.path(group_dir, "vi_halves_fast_summary.csv"), row.names = FALSE)
  print(summary_df)
}
