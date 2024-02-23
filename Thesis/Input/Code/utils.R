
# --------------- CLASS METHODS ----------------
setGeneric("flag", signature = "object",
           function(object) standardGeneric("flag"))


setMethod("flag", "SummarizedExperiment",
          function(object) rowData(object)$Flag)


setGeneric("flag<-", signature = "object",
           function(object, value) standardGeneric("flag<-"))


setMethod("flag<-", "SummarizedExperiment",
          function(object, value) {
            rowData(object)$Flag <- value
            if (validObject(object)) {
              return(object)
            }
          })

setGeneric("join_fData", signature = c("object", "dframe"),
           function(object, dframe) standardGeneric("join_fData"))

setMethod("join_fData", c("SummarizedExperiment", "data.frame"),
          function(object, dframe) {
            rowData(object) <- DataFrame(dplyr::left_join(as.data.frame(rowData(object)),
                                              dframe,
                                              by = "Feature_ID"))
            rownames(rowData(object)) <- rowData(object)$Feature_ID
            if (validObject(object)) {
              return(object)
            }
          })

# Class for creating a data frame with sample information and feature abundances as columns, with sample ID's on rows. For visualization functions.
setGeneric("combined_data", signature = "object",
           function(object) standardGeneric("combined_data"))

# Combined data method for SummarizedExperiment object
setMethod("combined_data", c(object = "SummarizedExperiment"),
          function(object) {
            cbind(as.data.frame(colData(object)), as.data.frame(t(assay(object))))
          })

# --------------- CLASS FUNCTIONS ----------------

drop_flagged <- function(object, all_features = FALSE) {
  if (!all_features) {
    object <- object[is.na(flag(object)), ]
  }
  object
}

drop_qcs <- function(object) {
  object <- object[, object$QC != "QC"]
  colData(object) <- droplevels(colData(object))
  object
}

mark_nas <- function(object, value, assay.type = NULL) {
  if (!is.null(assay.type)) {
    ex <- assays(object)[[assay.type]]
    ex[ex == value] <- NA
    assays(object)[[assay.type]] <- ex
  } else {
    ex <- assay(object)
    ex[ex == value] <- NA
    assay(object) <- ex
  }
  object
}

flag_quality <- function(object, assay.type = NULL,
                         condition = "(RSD_r < 0.2 & D_ratio_r < 0.4) | (RSD < 0.1 & RSD_r < 0.1 & D_ratio < 0.1)") {

  if (is.null(quality(object))) {
    object <- assess_quality(object, assay.type = assay.type)
  }
  good <- paste0("as.data.frame(rowData(object)) %>% dplyr::filter(", condition, ")") %>%
    parse(text = .) %>% eval()
  good <- good$Feature_ID

  idx <- is.na(flag(object)) & !rowData(object)$Feature_ID %in% good
  flag(object)[idx] <- "Low_quality"

  percentage <- scales::percent(sum(flag(object) == "Low_quality", na.rm = TRUE)/nrow(rowData(object)))

  object
}


# --------------- OTHER USEFUL FUNCTIONS ----------------

# Converter from MetaboSet to TreeSummarizedExperiment
library(TreeSummarizedExperiment)
makeTreeSummarizedExperimentFromMetaboSet <-
    function(object)
{
    assays <- exprs(object)
    colData <- DataFrame(pData(object))
    mcols(colData) <- DataFrame(varMetadata(object))
    rowData <- DataFrame(fData(object))
    mcols(rowData) <- DataFrame(varMetadata(featureData(object)))
    metadata <- SimpleList(
        experimentData = experimentData(object),
        annotation = annotation(object),
        protocolData = protocolData(object)
    )

    TreeSummarizedExperiment(
        assays = assays,
        colData = colData,
        rowData = rowData,
        metadata = metadata
    )
}

# Include p-values with sign of effect using instance, a rowData column specifying the sign of effect (direction_col) and the rowData column to be signed (x_col). For direction of effect with values around 1 (ratio, fold-change), use "type = 1", for direction of effect with values around 0 (difference in means), use "type = 0".
# sign_col(tse, direction_col = "diff_means", x_col = "pvalues", name= "sign_pvalues")
sign_col <- function(object, direction_col, x_col, type = 1, name){
  rowData(object)[, name] <- NA

  if(type == 0) {
    rowData(object)[rowData(object)[, direction_col] < 0, ][, name] <- rowData(object)[rowData(object)[, direction_col] < 0, ][, x_col] * -1

    rowData(object)[!rowData(object)[, direction_col] < 0, ][, name] <- rowData(object)[!rowData(object)[, direction_col] < 0, ][, x_col]
  }

  if(type == 1) {
    rowData(object)[rowData(object)[, direction_col] < 1, ][, name] <- rowData(object)[rowData(object)[, direction_col] < 1, ][, x_col] * -1

    rowData(object)[rowData(object)[, direction_col] >= 1, ][, name] <- rowData(object)[rowData(object)[, direction_col] >= 1, ][, x_col]
  }
  return(object)
}

# Make ranks from tiers based on other metric, for example fold change or p-value
rank_tiers <- function(object, tier_col, x_col) {
  rowData(object)$mul_ranks <- "NA"
  object <- object[order(rowData(object)[, x_col]), ]
  S <- grep("S", rowData(object)[, tier_col])
  A <- grep("A", rowData(object)[, tier_col])
  B <- grep("B", rowData(object)[, tier_col])
  C <- grep("C", rowData(object)[, tier_col])
  D <- grep("D", rowData(object)[, tier_col])
  E <- grep("E", rowData(object)[, tier_col])

  counter <- 1

  if (length(S) > 0) {
  rowData(object)[S, ]$mul_ranks <- counter:c(counter + c(length(S) -1))
  counter <- counter + length(S)
  }
  if (length(A) > 0) {
  rowData(object)[A, ]$mul_ranks <- counter:c(counter + c(length(A) -1))
  counter <- counter + length(A)
  }
  if (length(B) > 0) {
  rowData(object)[B, ]$mul_ranks <- counter:c(counter + c(length(B) -1))
  counter <- counter + length(B)
  }
  if (length(C) > 0) {
  rowData(object)[C, ]$mul_ranks <- counter:c(counter + c(length(C) -1))
  counter <- counter + length(C)
  }
  if (length(D) > 0) {
  rowData(object)[D, ]$mul_ranks <- counter:c(counter + c(length(D) -1))
  counter <- counter + length(D)
  }
  if (length(E) > 0) {
  rowData(object)[E, ]$mul_ranks <- 100000
  }
  rowData(object)$mul_ranks %<>% as.numeric()
  return(object)
}

# --------------- MINI HELPERS ----------------

prop_found <- function(x) {
  sum(!is.na(x)) / length(x)
}

finite_sd <- function(x) {
  sd(x[is.finite(x)], na.rm = TRUE)
}

finite_mean <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x[is.finite(x)], na.rm = TRUE)
}


finite_median <- function(x) {
  median(x[is.finite(x)], na.rm = TRUE)
}


finite_min <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  min(x[is.finite(x)], na.rm = TRUE)
}


finite_max <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  max(x[is.finite(x)], na.rm = TRUE)
}


finite_mad <- function(x) {
  mad(x[is.finite(x)], center = median(x[is.finite(x)], na.rm = TRUE), na.rm = TRUE)
}


finite_quantile <- function(x, ...) {
  unname(quantile(x[is.finite(x)], na.rm = TRUE, ...))
}

# Conditional pipe operator
"%||%" <- function(a, b) if (!is.null(a)) a else b

minus_log10 <- scales::trans_new("minus_log10",
                                 transform = function(x) {-log10(x)},
                                 inverse = function(x) {10^{-x}})
