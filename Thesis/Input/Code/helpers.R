fill_results <- function(results_df, features) {
  # Add NA rows for features where the test failed
  results_df <- results_df %>% dplyr::select(Feature_ID, dplyr::everything())
  missing_features <- setdiff(features, results_df$Feature_ID)
  fill_nas <- matrix(NA, nrow = length(missing_features), ncol = ncol(results_df) - 1) %>%
    as.data.frame()
  results_fill <- data.frame(Feature_ID = missing_features, fill_nas)
  rownames(results_fill) <- missing_features
  colnames(results_fill) <- colnames(results_df)
  results_df <- rbind(results_df, results_fill) %>% as.data.frame()
  rownames(results_df) <- results_df$Feature_ID
  # Set Feature ID to the original order
  results_df <- results_df[features, ]
  results_df
}

perform_test <- function(object, formula_char, result_fun, all_features, fdr = TRUE, packages = NULL) {

  data <- combined_data(object)
  features <- rownames(object)

  results_df <- foreach::foreach(i = seq_along(features), .combine = dplyr::bind_rows, .packages = packages) %dopar% {
    feature <- features[i]
    # Replace "Feature" with the current feature name
    tmp_formula <- gsub("Feature", feature, formula_char)
    # Run test
    result_row <- result_fun(feature = feature, formula = as.formula(tmp_formula), data = data)
    # In case Feature is used as predictor, make the column names match
    if (!is.null(result_row)){
      colnames(result_row) <- gsub(feature, "Feature", colnames(result_row))
    }
    result_row
  }
  # Check that results actually contain something
  # If the tests are run on parallel, the error messages from failing tests are not visible
  if (nrow(results_df) == 0) {
    stop("All the test failed, to see the individual error messages run the tests withot parallelization.",
         call. = FALSE)
  }
  # Rows full of NA for features where the test failed
  results_df <- fill_results(results_df, features)

  # FDR correction
  if (fdr) {
    if (all_features) {
      flags <- rep(NA_character_, nrow(results_df))
    } else {
      flags <- flag(object)
    }
    results_df <- adjust_p_values(results_df, flags)
  }
  results_df
}

density_plot <- function(data, x, fill, fill_scale = NULL, color_scale = NULL,
                         title = NULL, subtitle = NULL,
                         xlab = x, fill_lab = fill) {

  p <- ggplot(data, aes_string(x, fill = fill, color = fill)) +
    geom_density(alpha = 0.2) +
    fill_scale +
    labs(title = title, subtitle = subtitle, x = xlab, fill = fill_lab, color = NULL) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    color_scale

  p
}

quality <- function(object) {
  if (!all(c("RSD", "RSD_r", "D_ratio","D_ratio_r") %in% colnames(rowData(object)))) {
    return(NULL)
  }
  rowData(object)[c("Feature_ID", "RSD", "RSD_r", "D_ratio",
                  "D_ratio_r")]
}

erase_quality <- function(object) {
  if (!all(c("RSD", "RSD_r", "D_ratio","D_ratio_r") %in% colnames(rowData(object)))) {
    return(NULL)
  }
  rowData(object)[c("RSD", "RSD_r", "D_ratio", "D_ratio_r")] <- NULL
  object
}

assess_quality <- function(object, assay.type) {
  # Remove old quality metrics
  if (!is.null(quality(object))) {
    object <- erase_quality(object)
  }

  if (!is.null(assay.type)) {
    qc_data <- assays(object)[[assay.type]][, object$QC == "QC"]
    sample_data <- assays(object)[[assay.type]][, object$QC != "QC"]
  } else {
    qc_data <- assay(object)[, object$QC == "QC"]
    sample_data <- assay(object)[, object$QC != "QC"]
  }

  quality_metrics <- foreach::foreach(i = seq_len(nrow(sample_data)), .combine = rbind,
                                      .export = c("finite_sd", "finite_mad", "finite_mean", "finite_median")) %dopar% {
    data.frame(Feature_ID = rownames(sample_data)[i],
               RSD = finite_sd(qc_data[i, ]) / abs(finite_mean(qc_data[i, ])),
               RSD_r = finite_mad(qc_data[i, ]) / abs(finite_median(qc_data[i, ])),
               D_ratio = finite_sd(qc_data[i, ]) / finite_sd(sample_data[i, ]),
               D_ratio_r = finite_mad(qc_data[i, ]) / finite_mad(sample_data[i, ]),
               row.names = rownames(sample_data)[i], stringsAsFactors = FALSE)
  }

  object <- join_fData(object, quality_metrics)

  object
}

perform_lm <- function(object, formula_char, all_features = FALSE, ...) {

  lm_fun <- function(feature, formula, data) {
    # Try to fit the linear model
    fit <- NULL
    tryCatch({
      fit <- lm(formula, data = data, ...)
    }, error = function(e) cat(paste0(feature, ": ", e$message, "\n")))
    if(is.null(fit) | sum(!is.na(data[, feature])) < 2){
      result_row <- NULL
    } else {
      # Gather coefficients and CIs to one data frame row
      coefs <- summary(fit)$coefficients
      confints <- confint(fit, level = 0.95)
      coefs <- data.frame(Variable = rownames(coefs), coefs, stringsAsFactors = FALSE)
      confints <- data.frame(Variable = rownames(confints), confints, stringsAsFactors = FALSE)

      result_row <- dplyr::left_join(coefs,confints, by = "Variable") %>%
        dplyr::rename("Std_Error" = "Std..Error", "t_value" ="t.value",
                      "P" = "Pr...t..", "LCI95" = "X2.5..", "UCI95" = "X97.5..") %>%
        tidyr::gather("Metric", "Value", -Variable) %>%
        tidyr::unite("Column", Variable, Metric, sep="_") %>%
        tidyr::spread(Column, Value)
      # Add R2 statistics and feature ID
      result_row$R2 <- summary(fit)$r.squared
      result_row$Adj_R2 <- summary(fit)$adj.r.squared
      result_row$Feature_ID <- feature
    }
    result_row
  }

  results_df <- perform_test(object, formula_char, lm_fun, all_features)

  # Set a good column order
  variables <- gsub("_P$", "", colnames(results_df)[grep("P$", colnames(results_df))])
  statistics <- c("Estimate", "LCI95", "UCI95", "Std_Error", "t_value", "P", "P_FDR")
  col_order <- expand.grid(statistics, variables, stringsAsFactors = FALSE) %>%
    tidyr::unite("Column", Var2, Var1)
  col_order <- c("Feature_ID", col_order$Column, c("R2", "Adj_R2"))

  results_df[col_order]
}


adjust_p_values <- function(x, flags) {
  p_cols <- colnames(x)[grepl("_P$", colnames(x))]
  for (p_col in p_cols) {
    p_values <- x[, p_col, drop = TRUE]
    p_values[!is.na(flags)] <- NA
    x <- tibble::add_column(.data = x,
                            FDR = p.adjust(p_values, method = "BH"),
                            .after = p_col)
    p_idx <- which(colnames(x) == p_col)
    colnames(x)[p_idx + 1] <- paste0(p_col, "_FDR")
  }
  x
}
