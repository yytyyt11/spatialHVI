is_generic_sample_names <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(TRUE)
  }
  all(grepl("^V[0-9]+$", x))
}

strip_side_prefix <- function(x) {
  sub("^(Left|Right|left|right|L|R)[_.-]", "", x)
}

detect_halves_trait_order <- function(trait_names) {
  if (is.null(trait_names) || length(trait_names) == 0L || length(trait_names) %% 2L != 0L) {
    return("unknown")
  }

  is_left <- grepl("^(Left|left|L)[_.-]", trait_names)
  is_right <- grepl("^(Right|right|R)[_.-]", trait_names)
  stripped <- strip_side_prefix(trait_names)
  half <- length(trait_names) / 2L

  halves_ok <- all(is_left[seq_len(half)]) &&
    all(is_right[(half + 1L):length(trait_names)]) &&
    identical(stripped[seq_len(half)], stripped[(half + 1L):length(trait_names)])
  if (halves_ok) {
    return("halves")
  }

  odd_idx <- seq(1L, length(trait_names), by = 2L)
  even_idx <- seq(2L, length(trait_names), by = 2L)
  pairs_ok <- all(is_left[odd_idx]) &&
    all(is_right[even_idx]) &&
    identical(stripped[odd_idx], stripped[even_idx])
  if (pairs_ok) {
    return("pairs")
  }

  "unknown"
}

reorder_traits_to_halves <- function(y, trait_order = c("auto", "halves", "pairs")) {
  y <- validate_even_trait_matrix(y, arg = "y")
  trait_order <- match.arg(trait_order)
  trait_names <- rownames(y)
  detected_order <- detect_halves_trait_order(trait_names)
  applied_order <- trait_order

  if (trait_order == "auto") {
    applied_order <- if (detected_order == "pairs") "pairs" else "halves"
  }

  if (applied_order == "pairs") {
    reorder_idx <- c(seq(1L, nrow(y), by = 2L), seq(2L, nrow(y), by = 2L))
    y <- y[reorder_idx, , drop = FALSE]
    trait_names <- rownames(y)
  }

  list(
    y = y,
    trait_names = trait_names,
    detected_order = detected_order,
    applied_order = applied_order
  )
}

handle_missing_matrix <- function(mat, na_action = c("error", "zero"), arg = "mat") {
  na_action <- match.arg(na_action)
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"

  if (anyNA(mat)) {
    na_count <- sum(is.na(mat))
    if (na_action == "error") {
      stop(sprintf("`%s` contains %d missing values.", arg, na_count), call. = FALSE)
    }
    mat[is.na(mat)] <- 0
  }

  mat[!is.finite(mat)] <- 0
  mat
}

read_brain_measure_input <- function(
  brain_csv,
  brain_row_names = TRUE,
  phenotype_layout = c("auto", "traits_in_rows", "samples_in_rows"),
  phenotype_id_col = NULL,
  phenotype_sample_ids = NULL,
  trait_order = c("auto", "halves", "pairs"),
  na_action = c("error", "zero")
) {
  phenotype_layout <- match.arg(phenotype_layout)
  trait_order <- match.arg(trait_order)
  na_action <- match.arg(na_action)

  if (!file.exists(brain_csv)) {
    stop(sprintf("File not found: %s", brain_csv), call. = FALSE)
  }

  df <- if (isTRUE(brain_row_names)) {
    utils::read.csv(
      brain_csv,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      row.names = 1
    )
  } else {
    utils::read.csv(
      brain_csv,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (phenotype_layout == "auto") {
    row_order <- detect_halves_trait_order(rownames(df))
    col_order <- detect_halves_trait_order(colnames(df))
    phenotype_layout <- if (row_order != "unknown") {
      "traits_in_rows"
    } else if (col_order != "unknown") {
      "samples_in_rows"
    } else if (isTRUE(brain_row_names)) {
      "traits_in_rows"
    } else {
      "samples_in_rows"
    }
  }

  if (phenotype_layout == "traits_in_rows") {
    non_numeric <- names(df)[!vapply(df, is.numeric, logical(1))]
    if (length(non_numeric) > 0L) {
      stop(
        sprintf(
          "All phenotype columns must be numeric when `phenotype_layout = \"traits_in_rows\"`; non-numeric columns: %s",
          paste(non_numeric, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    y <- handle_missing_matrix(df, na_action = na_action, arg = "brain_csv")
    sample_ids <- if (!is.null(phenotype_sample_ids)) {
      as.character(phenotype_sample_ids)
    } else if (is_generic_sample_names(colnames(y))) {
      NULL
    } else {
      colnames(y)
    }
    if (!is.null(sample_ids) && length(sample_ids) != ncol(y)) {
      stop("`phenotype_sample_ids` must have length equal to the number of samples.", call. = FALSE)
    }

    reordered <- reorder_traits_to_halves(y, trait_order = trait_order)
    return(
      list(
        y = reordered$y,
        raw_y = y,
        trait_names = rownames(reordered$y),
        sample_ids = sample_ids,
        phenotype_layout = phenotype_layout,
        trait_order_detected = reordered$detected_order,
        trait_order_applied = reordered$applied_order
      )
    )
  }

  id_values <- NULL
  if (!is.null(phenotype_id_col)) {
    if (length(phenotype_id_col) == 1L && is.character(phenotype_id_col)) {
      if (!phenotype_id_col %in% names(df)) {
        stop(sprintf("`phenotype_id_col` `%s` was not found.", phenotype_id_col), call. = FALSE)
      }
      id_values <- as.character(df[[phenotype_id_col]])
      df <- df[, setdiff(names(df), phenotype_id_col), drop = FALSE]
    } else {
      stop("`phenotype_id_col` must be a single column name.", call. = FALSE)
    }
  } else if (!is.null(rownames(df)) && !all(rownames(df) %in% as.character(seq_len(nrow(df))))) {
    id_values <- rownames(df)
  }

  non_numeric <- names(df)[!vapply(df, is.numeric, logical(1))]
  if (length(non_numeric) > 0L) {
    stop(
      sprintf(
        "All phenotype trait columns must be numeric after removing the ID column; non-numeric columns: %s",
        paste(non_numeric, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  y <- handle_missing_matrix(t(as.matrix(df)), na_action = na_action, arg = "brain_csv")
  reordered <- reorder_traits_to_halves(y, trait_order = trait_order)

  list(
    y = reordered$y,
    raw_y = y,
    trait_names = rownames(reordered$y),
    sample_ids = id_values,
    phenotype_layout = phenotype_layout,
    trait_order_detected = reordered$detected_order,
    trait_order_applied = reordered$applied_order
  )
}

auto_detect_unique_id_col <- function(df) {
  candidate_cols <- names(df)[vapply(df, function(col) !is.numeric(col), logical(1))]
  if (length(candidate_cols) == 0L) {
    return(NULL)
  }

  unique_cols <- candidate_cols[vapply(
    df[candidate_cols],
    function(col) !anyNA(col) && length(unique(col)) == length(col),
    logical(1)
  )]
  if (length(unique_cols) == 0L) {
    return(NULL)
  }
  if ("PTID" %in% unique_cols) {
    return("PTID")
  }
  unique_cols[1]
}

resolve_column_selection <- function(x, selection, arg = "selection") {
  if (is.null(selection)) {
    return(NULL)
  }
  if (is.character(selection)) {
    missing <- setdiff(selection, colnames(x))
    if (length(missing) > 0L) {
      stop(sprintf("Unknown column names in `%s`: %s", arg, paste(missing, collapse = ", ")), call. = FALSE)
    }
    return(selection)
  }
  if (is.numeric(selection)) {
    return(colnames(x)[selection])
  }
  stop(sprintf("`%s` must be NULL, column names, or column indices.", arg), call. = FALSE)
}

load_rdata_env <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path), call. = FALSE)
  }

  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  list(env = env, objects = loaded)
}

detect_rdata_object_name <- function(env, object_names, kind = c("x", "k", "annotation")) {
  kind <- match.arg(kind)
  candidates <- character(0)

  for (nm in object_names) {
    obj <- get(nm, envir = env, inherits = FALSE)
    if (kind == "annotation") {
      if (is.data.frame(obj) && "RSID" %in% names(obj)) {
        candidates <- c(candidates, nm)
      }
      next
    }

    if (!(is.matrix(obj) || is.data.frame(obj))) {
      next
    }

    obj_df <- as.data.frame(obj, stringsAsFactors = FALSE)
    numeric_cols <- names(obj_df)[vapply(obj_df, is.numeric, logical(1))]
    if (kind == "x") {
      id_col <- auto_detect_unique_id_col(obj_df)
      snp_like <- grepl("^rs", names(obj_df))
      if (!is.null(id_col) || any(snp_like) || length(numeric_cols) > 10L) {
        candidates <- c(candidates, nm)
      }
    } else if (kind == "k") {
      if (nrow(obj_df) == ncol(obj_df) && length(numeric_cols) == ncol(obj_df)) {
        mat <- as.matrix(obj_df)
        if (isTRUE(all.equal(mat, t(mat)))) {
          candidates <- c(candidates, nm)
        }
      }
    }
  }

  if (length(candidates) == 0L) {
    return(NULL)
  }
  candidates[1]
}

prepare_realdata_x <- function(
  x,
  x_ids = NULL,
  x_id_col = NULL,
  x_snp_cols = NULL,
  x_snp_pattern = "^rs",
  annotation = NULL,
  annotation_rsid_col = "RSID"
) {
  if (is.data.frame(x)) {
    x_df <- x
  } else if (is.matrix(x)) {
    x_df <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
    names(x_df) <- colnames(x)
  } else {
    stop("`X` must be a matrix or data frame.", call. = FALSE)
  }

  if (is.null(x_id_col)) {
    x_id_col <- auto_detect_unique_id_col(x_df)
  }
  if (!is.null(x_id_col)) {
    if (!x_id_col %in% names(x_df)) {
      stop(sprintf("`x_id_col` `%s` was not found in `X`.", x_id_col), call. = FALSE)
    }
    sample_ids <- as.character(x_df[[x_id_col]])
    x_df <- x_df[, setdiff(names(x_df), x_id_col), drop = FALSE]
  } else if (!is.null(x_ids)) {
    sample_ids <- as.character(x_ids)
  } else if (!is.null(rownames(x_df)) && !all(rownames(x_df) %in% as.character(seq_len(nrow(x_df))))) {
    sample_ids <- rownames(x_df)
  } else {
    sample_ids <- NULL
  }

  chosen_cols <- resolve_column_selection(x_df, x_snp_cols, arg = "x_snp_cols")
  if (is.null(chosen_cols) && !is.null(annotation) && annotation_rsid_col %in% names(annotation)) {
    annotation_ids <- as.character(annotation[[annotation_rsid_col]])
    chosen_cols <- names(x_df)[names(x_df) %in% annotation_ids]
  }
  if (is.null(chosen_cols)) {
    pattern_hits <- grepl(x_snp_pattern, names(x_df))
    if (any(pattern_hits)) {
      chosen_cols <- names(x_df)[pattern_hits]
    }
  }
  if (is.null(chosen_cols)) {
    chosen_cols <- names(x_df)[vapply(x_df, is.numeric, logical(1))]
  }
  if (length(chosen_cols) == 0L) {
    stop("No SNP columns were identified in `X`.", call. = FALSE)
  }

  non_numeric <- chosen_cols[!vapply(x_df[chosen_cols], is.numeric, logical(1))]
  if (length(non_numeric) > 0L) {
    stop(
      sprintf("Selected SNP columns must be numeric; non-numeric columns: %s", paste(non_numeric, collapse = ", ")),
      call. = FALSE
    )
  }

  extra_cols <- setdiff(names(x_df), chosen_cols)
  extra_numeric_cols <- extra_cols[vapply(x_df[extra_cols], is.numeric, logical(1))]
  extra_numeric <- x_df[, extra_numeric_cols, drop = FALSE]
  x_mat <- handle_missing_matrix(as.matrix(x_df[chosen_cols]), na_action = "error", arg = "X")

  if (!is.null(sample_ids) && length(sample_ids) != nrow(x_mat)) {
    stop("Sample IDs extracted from `X` do not match the number of rows.", call. = FALSE)
  }

  list(
    X = x_mat,
    sample_ids = sample_ids,
    snp_names = chosen_cols,
    shared_numeric = extra_numeric,
    id_col = x_id_col
  )
}

prepare_realdata_k <- function(K, k_ids = NULL) {
  if (is.data.frame(K)) {
    k_df <- K
  } else if (is.matrix(K)) {
    k_df <- as.data.frame(K, stringsAsFactors = FALSE, check.names = FALSE)
    names(k_df) <- colnames(K)
  } else {
    stop("`K` must be a matrix or data frame.", call. = FALSE)
  }

  id_col <- auto_detect_unique_id_col(k_df)
  if (!is.null(id_col)) {
    sample_ids <- as.character(k_df[[id_col]])
    k_df <- k_df[, setdiff(names(k_df), id_col), drop = FALSE]
  } else if (!is.null(k_ids)) {
    sample_ids <- as.character(k_ids)
  } else if (!is.null(rownames(k_df)) && !all(rownames(k_df) %in% as.character(seq_len(nrow(k_df))))) {
    sample_ids <- rownames(k_df)
  } else {
    sample_ids <- NULL
  }

  non_numeric <- names(k_df)[!vapply(k_df, is.numeric, logical(1))]
  if (length(non_numeric) > 0L) {
    stop(
      sprintf("All kinship columns must be numeric; non-numeric columns: %s", paste(non_numeric, collapse = ", ")),
      call. = FALSE
    )
  }

  k_mat <- handle_missing_matrix(as.matrix(k_df), na_action = "error", arg = "K")
  k_mat <- validate_square_matrix(k_mat, arg = "K")

  if (!is.null(sample_ids) && length(sample_ids) != nrow(k_mat)) {
    sample_ids <- NULL
  }
  if (is.null(sample_ids) && !is.null(rownames(k_mat)) && !is.null(colnames(k_mat)) && identical(rownames(k_mat), colnames(k_mat))) {
    sample_ids <- rownames(k_mat)
  }

  list(K = k_mat, sample_ids = sample_ids)
}

validate_shared_trait_alignment <- function(y, trait_names, shared_numeric, tolerance = 1e-8) {
  if (is.null(shared_numeric) || ncol(shared_numeric) == 0L || is.null(trait_names)) {
    return(list(ok = NA, matched_traits = character(0), max_abs_diff = numeric(0)))
  }

  shared_traits <- intersect(trait_names, names(shared_numeric))
  if (length(shared_traits) == 0L) {
    return(list(ok = NA, matched_traits = character(0), max_abs_diff = numeric(0)))
  }

  if (ncol(y) != nrow(shared_numeric)) {
    return(list(ok = FALSE, matched_traits = shared_traits, max_abs_diff = Inf))
  }

  diffs <- vapply(
    shared_traits,
    function(nm) max(abs(as.numeric(y[nm, ]) - as.numeric(shared_numeric[[nm]]))),
    numeric(1)
  )
  list(ok = all(diffs <= tolerance), matched_traits = shared_traits, max_abs_diff = diffs)
}

align_realdata_samples <- function(
  y,
  phenotype_sample_ids = NULL,
  x_mat = NULL,
  x_ids = NULL,
  k_mat = NULL,
  k_ids = NULL,
  allow_row_order_alignment = TRUE,
  shared_validation = NULL
) {
  phenotype_sample_ids <- if (is.null(phenotype_sample_ids)) NULL else as.character(phenotype_sample_ids)
  x_ids <- if (is.null(x_ids)) NULL else as.character(x_ids)
  k_ids <- if (is.null(k_ids)) NULL else as.character(k_ids)

  if (!is.null(x_mat)) {
    if (!is.null(phenotype_sample_ids) && !is.null(x_ids)) {
      common_ids <- intersect(phenotype_sample_ids, x_ids)
      if (length(common_ids) < 2L) {
        stop("Phenotype sample IDs and `X` sample IDs do not overlap enough to align.", call. = FALSE)
      }
      y <- y[, match(common_ids, phenotype_sample_ids), drop = FALSE]
      x_mat <- x_mat[match(common_ids, x_ids), , drop = FALSE]
      return(list(y = y, x_mat = x_mat, sample_ids = common_ids, method = "sample_ids"))
    }

    if (is.list(shared_validation) && identical(shared_validation$ok, TRUE)) {
      if (ncol(y) != nrow(x_mat)) {
        stop("Phenotype columns and rows of `X` do not have the same sample count.", call. = FALSE)
      }
      return(list(y = y, x_mat = x_mat, sample_ids = x_ids, method = "shared_traits"))
    }

    if (is.list(shared_validation) && identical(shared_validation$ok, FALSE)) {
      stop("Shared phenotype columns in `X` do not match the sample order in `brain_csv`.", call. = FALSE)
    }

    if (!isTRUE(allow_row_order_alignment)) {
      stop(
        paste(
          "No sample IDs were available to align `brain_csv` and `X`.",
          "Set `allow_row_order_alignment = TRUE` only if their order already matches."
        ),
        call. = FALSE
      )
    }
    if (ncol(y) != nrow(x_mat)) {
      stop("Phenotype columns and rows of `X` do not have the same sample count.", call. = FALSE)
    }
    warning(
      "Aligned phenotype and genotype data by row order because no usable sample IDs were available.",
      call. = FALSE
    )
    return(list(y = y, x_mat = x_mat, sample_ids = x_ids, method = "row_order"))
  }

  if (!is.null(k_mat)) {
    if (!is.null(phenotype_sample_ids) && !is.null(k_ids)) {
      common_ids <- intersect(phenotype_sample_ids, k_ids)
      if (length(common_ids) < 2L) {
        stop("Phenotype sample IDs and `K` sample IDs do not overlap enough to align.", call. = FALSE)
      }
      y <- y[, match(common_ids, phenotype_sample_ids), drop = FALSE]
      idx <- match(common_ids, k_ids)
      k_mat <- k_mat[idx, idx, drop = FALSE]
      return(list(y = y, k_mat = k_mat, sample_ids = common_ids, method = "sample_ids"))
    }

    if (!isTRUE(allow_row_order_alignment)) {
      stop(
        paste(
          "No sample IDs were available to align `brain_csv` and `K`.",
          "Set `allow_row_order_alignment = TRUE` only if their order already matches."
        ),
        call. = FALSE
      )
    }
    if (ncol(y) != nrow(k_mat)) {
      stop("Phenotype columns and dimensions of `K` do not have the same sample count.", call. = FALSE)
    }
    warning(
      "Aligned phenotype and kinship data by row order because no usable sample IDs were available.",
      call. = FALSE
    )
    return(list(y = y, k_mat = k_mat, sample_ids = k_ids, method = "row_order"))
  }

  stop("Either `x_mat` or `k_mat` must be supplied for alignment.", call. = FALSE)
}
