if (!file.exists("BrainMeasure_data.csv") || !file.exists("FreeSurfer_Data.RData")) {
  stop("BrainMeasure_data.csv and FreeSurfer_Data.RData must exist in the package root.")
}

devtools::load_all(quiet = TRUE)

brain <- read.csv("BrainMeasure_data.csv", row.names = 1, check.names = FALSE)
brain_small <- brain[, 1:80, drop = FALSE]
brain_small_path <- file.path(tempdir(), "brain_small.csv")
write.csv(brain_small, brain_small_path, row.names = TRUE)

fit <- run_realdata_halves_from_files(
  brain_csv = brain_small_path,
  rdata_path = "FreeSurfer_Data.RData",
  use_k_from = "X",
  x_object = "FreeSurfer_list_Data",
  annotation_object = "SNP_gene_member_reduced",
  x_id_col = "PTID",
  x_snp_cols = 1:50,
  genotype_rows = 1:80,
  brain_row_names = TRUE,
  phenotype_layout = "traits_in_rows",
  max_iter = 1,
  tol = 0,
  verbose = FALSE,
  prefix = NULL,
  make_plot = FALSE,
  show_summary = FALSE
)

cat("realdata smoke test summary\n")
cat("alignment_method:", fit$alignment_method, "\n")
cat("sample_count:", length(fit$sample_ids), "\n")
cat("trait_count:", length(fit$trait_names), "\n")
cat("shared_validation_ok:", fit$shared_trait_validation$ok, "\n")
cat("E_rho:", fit$res$E_rho, "\n")
cat("H2:", fit$h$H2, "\n")
