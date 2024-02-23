# Tse instances for testing visualizations
# tse <- readRDS("post_drift_tse.RDS")
# intestine_tse <- readRDS("dummy_intestine_tse.RDS")


# ----------------- QC VISUALIZATIONS -----------------

# injection_lm_plot(tse)
injection_lm_plot <- function(object, assay.type = NULL, title = NULL,
                              subtitle = NULL) {

  if (!is.null(assay.type)) {
    assays(object)[!names(assays(object)) %in% assay.type] <- NULL
  }
  lm_all <- perform_lm(object, "Feature ~ Injection_order")
  lm_sample <- perform_lm(object[, object$QC != "QC"], "Feature ~ Injection_order")
  lm_qc <- perform_lm(object[, object$QC == "QC"], "Feature ~ Injection_order")

  p_values <- list("All samples" = lm_all$Injection_order_P,
                   "Biological samples" = lm_sample$Injection_order_P,
                   "QC samples" = lm_qc$Injection_order_P)
  p_histogram_plot(p_values, title=title, subtitle=subtitle)
}


# sample_boxplot_plot(tse, order_by="Group", fill_by="Group")
sample_boxplot_plot <- function(object, assay.type = NULL, order_by, fill_by,
                                title = NULL, subtitle = NULL,
                                fill_scale = ggplot2::scale_fill_brewer(palette = "Set1"), zoom_boxplot = TRUE) {
  if (!is.null(assay.type)) {
    assays(object)[!names(assays(object)) %in% assay.type] <- NULL
  }

  data <- combined_data(object)

  if (length(order_by) == 1) {
    data$order_by <- data[, order_by]
  } else {
    data <- tidyr::unite(data, "order_by", order_by, remove = FALSE)
  }

  if (length(fill_by) == 1) {
    data$fill_by <- data[, fill_by]
  } else {
    data <- tidyr::unite(data, "fill_by", fill_by, remove = FALSE)
  }

  data <- data %>%
    dplyr::arrange(order_by)

  data$Sample_ID <- factor(data$Sample_ID, levels = data$Sample_ID)

  data <- tidyr::gather(data, "Variable", "Value", rownames(rowData(object)))

  p <- ggplot(data, aes(x = Sample_ID, y = Value, fill = fill_by))

  if(zoom_boxplot){
    ylimits <- data %>%
      dplyr::group_by(Sample_ID) %>%
      dplyr::summarise(low = boxplot.stats(Value)$stats[1],
                       high = boxplot.stats(Value)$stats[5])

    ylimits <- c(0, max(ylimits$high))
    p <- p +
      geom_boxplot(outlier.shape = NA, na.rm=TRUE) +
      coord_cartesian(ylim = ylimits)
    subtitle <- paste(subtitle, "(zoomed in boxplot: outliers out of view)", sep=" ")
  } else {
    p <- p +
      geom_boxplot(na.rm=TRUE)
  }

  p <- p +
    labs(x = paste(order_by, collapse = "_"),
         y = "Abundance of metabolites",
         fill = paste(fill_by, collapse = "_"),
         title = title, subtitle = subtitle) +
    theme_bw() +
    theme(plot.title = element_text(face="bold"),
          axis.text.x = element_text(angle=90, vjust=0.3)) +
    fill_scale

  p
}

# dist_density_plot(tse)
dist_density_plot <- function(object, assay.type = NULL, all_features = FALSE,
                              dist_method = "euclidean", center = TRUE,
                              scale = "uv",
                              color_scale = ggplot2::scale_color_brewer(palette = "Set1"),
                              fill_scale = ggplot2::scale_fill_brewer(palette = "Set1"),
                              title = paste("Density plot of", dist_method, "distances between samples"),
                              subtitle = NULL) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
      stop("Package \"pcaMethods\" needed for this function to work. Please install it.",
           call. = FALSE)
  }

  if (!is.null(assay.type)) {
    assays(object)[!names(assays(object)) %in% assay.type] <- NULL
  }

  assay(object) <- pcaMethods::prep(assay(object), center = center, scale = scale)

  qc_data <- t(assay(object)[, object$QC == "QC"])
  sample_data <- t(assay(object)[, object$QC != "QC"])

  qc_dist <- dist(qc_data, method = dist_method) %>% as.numeric()
  sample_dist <- dist(sample_data, method = dist_method) %>% as.numeric()
  qc <- rep(c("QC", "Sample"), times = c(length(qc_dist), length(sample_dist)))
  qc <- rep(c("QC", "Sample"), times = c(length(qc_dist), length(sample_dist)))
  distances <- data.frame(dist = c(qc_dist, sample_dist), qc = qc)

  density_plot(distances, x = "dist", fill = "qc", fill_scale = fill_scale,
               color_scale = color_scale,
               xlab = "Distance", fill_lab = NULL,
               title = title, subtitle = subtitle)
}

# dendrogram_plot(tse, color="Group")
dendrogram_plot <- function(object, assay.type = NULL, color,
                            dist_method = "euclidean",
                            clust_method = "ward.D2", center = TRUE,
                            scale = "uv",
                            title = "Dendrogram of hierarchical clustering",
                            subtitle = NULL,
                            color_scale = ggplot2::scale_color_brewer(palette = "Set1")) {

  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
      stop("Package \"pcaMethods\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  if (!requireNamespace("ggdendro", quietly = TRUE)) {
      stop("Package \"ggdendro\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  color <- color %||% NULL

  subtitle <- subtitle %||% paste("Distance method:", dist_method, "Clustering method:", clust_method)

  if (!is.null(assay.type)) {
    assays(object)[!names(assays(object)) %in% assay.type] <- NULL
  }

  assay(object) <- pcaMethods::prep(assay(object), center = center, scale = scale)

  d_data <- dist(t(assay(object)), method = dist_method) %>%
    hclust(method = clust_method) %>%
    as.dendrogram() %>%
    ggdendro::dendro_data()

  labels <- ggdendro::label(d_data) %>%
    dplyr::mutate(label = as.character(label)) %>%
    dplyr::left_join(as.data.frame(colData(object))[c("Sample_ID", color)], by = c("label" = "Sample_ID"))
  labels[, color] <- as.factor(labels[, color])

  p <- ggplot(ggdendro::segment(d_data)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data = labels, aes_string(x = "x", y = "y", label = "label", color = color), angle = 90, hjust = 1) +
    ggdendro::theme_dendro() +
    color_scale +
    labs(title = title, subtitle = subtitle)

  p
}

# sample_heatmap_plot(tse, group_col = "Group")
sample_heatmap_plot <- function(object, assay.type = NULL, group_col,
                                group = colData(object)[, group_col],
                                dist_method = "euclidean",
                                clust_method = "ward.D2", center = TRUE,
                                scale = "uv", group_bar = TRUE,
                                title = "Heatmap of distances between samples",
                                subtitle = NULL,
                                fill_scale_con = ggplot2::scale_fill_viridis_c(),
                                fill_scale_dis = ggplot2::scale_fill_brewer(palette = "Set1")) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
      stop("Package \"pcaMethods\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  if (!requireNamespace("cowplot", quietly = TRUE)) {
      stop("Package \"cowplot\" needed for this function to work. Please install it.",
           call. = FALSE)
  }

  subtitle <- subtitle %||% paste("Distance method:", dist_method, "Clustering method:", clust_method)

  if (!is.null(assay.type)) {
    assays(object)[!names(assays(object)) %in% assay.type] <- NULL
  }

  assay(object) <- pcaMethods::prep(assay(object), center = center, scale = scale)

  distances <- dist(t(as.data.frame(assay(object))), method = dist_method)

  hc <- hclust(distances, method = clust_method)
  hc_order <- hc$labels[hc$order]

  distances_df <- as.matrix(distances) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("X") %>%
    tidyr::gather("Y", "Distance", -X)
  distances_df$X <- factor(distances_df$X, levels = hc_order, ordered = TRUE)
  distances_df$Y <- factor(distances_df$Y, levels = rev(hc_order), ordered = TRUE)

  p <- ggplot(distances_df, aes(X, Y, fill = Distance)) +
    geom_tile(color = NA) +
    fill_scale_con +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.05)) +
    coord_fixed()

  pheno_data <- as.data.frame(colData(object))
  pheno_data$Sample_ID <- factor(pheno_data$Sample_ID, levels = hc_order, ordered = TRUE)

  gb <- ggplot(pheno_data, aes(x = Sample_ID, y = 1, fill = group)) +
      geom_tile(color = "white") +
      theme_void() +
      fill_scale_dis

  p <- cowplot::plot_grid(p, gb, ncol = 1, align = "v", rel_heights = c(10/11,1/11))

  p
}

# --------------- FEATURE-WISE RESULTS VISUALIZATIONS ----------------

# beeswarm_plot(intestine_tse[1], assay.type = "aggregated", group_col="Group", color = "Group", add_boxplots = TRUE)
beeswarm_plot <- function(object, assay.type = NULL, group_col, color,
                          add_boxplots = FALSE,
                          color_scale=ggplot2::scale_color_brewer(palette = "Set1")) {

  if (!is.null(assay.type)) {
    assays(object)[!names(assays(object)) %in% assay.type] <- NULL
  }

  data <- combined_data(object)

  p <- ggplot(data, aes_string(x = group_col, y = rownames(object), color = color)) +
  theme_bw() +
  labs(y="Abundance") +
  ggtitle(rownames(object))
  if (add_boxplots) {
    p <- p +
      geom_boxplot(position = position_dodge(0.6), width = 0.5, lwd = .3) +
      stat_boxplot(geom ='errorbar', width = 0.5, lwd = .3)
  }
  p <- p +
    ggbeeswarm::geom_beeswarm() +
    color_scale
  p
}

# group_boxplot_plot(intestine_tse[1], assay.type = "aggregated",  group_col="Group", color="Group")
group_boxplot_plot <- function(object, assay.type = NULL, group_col, color,
                               color_scale = ggplot2::scale_color_brewer(palette = "Set1")) {

  if (!is.null(assay.type)) {
    assays(object)[!names(assays(object)) %in% assay.type] <- NULL
  }

  data <- combined_data(object)

  p <- ggplot(data, aes_string(x = group_col, y = rownames(object), color = color)) +
    labs(y="Abundance") +
    theme_bw() +
    ggtitle(rownames(object)) +
    geom_boxplot(position = position_dodge(0.6), width = 0.5) +
    stat_summary(fun.data = mean_se,
                 geom = "point",
                 shape = 18,
                 size = 3,
                 position = position_dodge(0.6)
    ) +
    color_scale
  p
}

# subject_line_plot(intestine_tse[1], assay.type="aggregated", time_col="Time", id_col="Subject_ID", color="Group")
subject_line_plot <- function(object, assay.type = NULL, time_col, id_col,
                              color=NA, facet=NULL, color_scale=ggplot2::scale_color_brewer(palette = "Set1")) {

  if (!is.null(assay.type)) {
    assays(object)[!names(assays(object)) %in% assay.type] <- NULL
  }

  data <- combined_data(object)

  p <- ggplot(data, aes_string(x = time_col, y = rownames(object))) +
    labs(y="Abundance") +
    theme_bw() +
    ggtitle(rownames(object))

  if (is.na(color)) {
    p <- p +
      geom_line(aes_string(group = id_col),
                color = "grey20",
                alpha = 0.35,
                size = 0.3
      ) +
      stat_summary(aes(group = 1),
                    fun.data = "mean_se",
                    geom = "line",
                    size = 1.2,
                    color = "red"
      )
  } else {
    p <- p +
      geom_line(aes_string(group = id_col, color = color),
                alpha = 0.35,
                size = 0.3
      ) +
      stat_summary(aes_string(group = color, color = color),
                    fun.data = "mean_se",
                    geom = "line",
                    size = 1.2
      ) +
      color_scale
  }

  if (!is.null(facet)) {
    p <- p + facet_wrap(facets = facet)
  }

  if (class(data[, time_col]) == "factor") {
    p <- p +
      scale_x_discrete(expand = c(0.05,0.05))
  }
  p
}

# group_line_plot(intestine_tse[1], assay.type="aggregated", time_col="Time", group_col="Group")
group_line_plot <- function(object, assay.type = NULL, time_col, group_col,
                            fun.data = "mean_cl_boot",
                            position_dodge_amount = 0.2, color_scale=ggplot2::scale_color_brewer(palette = "Set1"),
                            fun = NULL, fun.min = NULL, fun.max = NULL) {

  if (!is.null(assay.type)) {
    assays(object)[!names(assays(object)) %in% assay.type] <- NULL
  }

  data <- combined_data(object)

  p <- ggplot(data,
              aes_string(x = time_col, y = rownames(object), group =group_col, color = group_col)
  ) +
    theme_bw() +
    ggtitle(rownames(object)) +
    stat_summary(fun.data = fun.data,
                  geom = "errorbar", width = 0.5,
                  fun = fun,
                  fun.min = fun.min,
                  fun.max = fun.max,
                  position = position_dodge(position_dodge_amount)
    ) +
    stat_summary(fun.data = fun.data,
                  geom = "point",
                  fun = fun,
                  fun.min = fun.min,
                  fun.max = fun.max,
                  position = position_dodge(position_dodge_amount),
                  size = 4
    ) +
    stat_summary(fun.data = fun.data,
                  geom = "line",
                  position = position_dodge(position_dodge_amount), size = 0.5,
                  fun = fun,
                  fun.min = fun.min,
                  fun.max = fun.max
    ) +
    color_scale
    p
}

# --------------- COMPREHENSIVE RESULTS VISUALIZATIONS ----------------

# p_histogram_plot(p_values=list(rowData(intestine_tse)$pvalues))
p_histogram_plot <- function(p_values, hline = TRUE, combine = TRUE,
                             x_label = "p-value", title = NULL,
                             subtitle = NULL) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
      stop("Package \"cowplot\" needed for this function to work. Please install it.",
           call. = FALSE)
  }
  breaks <- seq(0, 1, by = 0.05)

  plots <- list()
  for (i in seq_along(p_values)) {

    p <- ggplot(data.frame(P = p_values[[i]]), aes(P)) +
      geom_histogram(breaks = breaks, col = "grey50", fill = "grey80", size = 1) +
      labs(x = x_label, y = "Frequency") +
      ggtitle(names(p_values)[i]) +
      theme_minimal() +
      theme(plot.title = element_text(face="bold", hjust=0.5))

    if (hline) {
      finite_count <- sum(is.finite(p_values[[i]]))
      h_line <- finite_count/(length(breaks)-1)
      p <- p +
        geom_hline(yintercept = h_line, color="red", linetype = "dashed", size = 1)
    }

    plots <- c(plots, list(p))
  }

  title_gg <- ggplot() +
              labs(title = title, subtitle = subtitle)
  if (combine) {
    grid <- cowplot::plot_grid(plotlist = plots, ncol = 1)
    if (!is.null(title)) {
      p <- cowplot::plot_grid(title_gg, grid, ncol = 1, rel_heights = c(0.05, 1))
    } else {
      p <- cowplot::plot_grid(plotlist = plots, ncol = 1)
    }
  } else {
    if (!is.null(title)) {
      p <- cowplot::plot_grid(title_gg, plotlist=plots, ncol = 1, rel_heights = c(0.05, 1))
    }
  }
  return(p)
}


# volcano_plot(intestine_tse, effect_col="FC", p_col="pvalues", q_col="pvalues", color="pvalues", log2_x=TRUE)
setGeneric("volcano_plot", signature = "object",
           function(object, effect_col, p_col, q_col = NULL, color = NULL,
                    p_breaks = c(0.05, 0.01, 0.001, 1e-4), fdr_limit = 0.05,
                    log2_x = FALSE, center_x_axis = TRUE, x_lim = NULL, label = NULL, label_limit = 0.05,
                    color_scale = ggplot2::scale_color_viridis_c(),
                    title = "Volcano plot", subtitle = NULL, ...) standardGeneric("volcano_plot"))

setMethod("volcano_plot", c(object = "SummarizedExperiment"),
          function(object, effect_col, p_col, q_col = NULL, color = NULL,
                   p_breaks = c(0.05, 0.01, 0.001, 1e-4), fdr_limit = 0.05,
                   log2_x = FALSE, center_x_axis = TRUE, x_lim = NULL, label = NULL, label_limit = 0.05,
                   color_scale = ggplot2::scale_color_viridis_c(),
                   title = "Volcano plot", subtitle = NULL, ...) {
            volcano_plotter(as.data.frame(rowData(object)), effect_col, p_col,
                            q_col, color, p_breaks, fdr_limit, log2_x, center_x_axis, x_lim, label, label_limit, color_scale, title, subtitle, ...)
          })

volcano_plotter <- function(data, effect_col, p_col, q_col, color, p_breaks,
                            fdr_limit,log2_x, center_x_axis, x_lim, label, label_limit, color_scale, title, subtitle, ...) {

  if (center_x_axis & !is.null(x_lim)) {
    warning("Manually setting x-axis limits overrides x-axis centering")
    center_x_axis <- FALSE
  }
  if (min(data[, p_col]) > max(p_breaks)) {
    warning("All the p-values are larger than the p-value breaks supplied. Consider using larger p_breaks for plotting")
  }

  pl <- ggplot(data, aes_string(x = effect_col, y = p_col, color = color)) +
    geom_point(...) +
    color_scale +
    theme_bw() +
    labs(title = title, subtitle = subtitle, y="p-value") +
    theme(panel.grid.minor.y = element_blank(),
          axis.ticks.y = element_blank())

  if(!is.null(q_col)) {

    if (any(data[, q_col] < fdr_limit)) {
      q_limit <- max(data[data[, q_col] < fdr_limit, p_col], na.rm = TRUE)
      pl <- pl +
        geom_hline(yintercept = q_limit, linetype = "dashed") +
        scale_y_continuous(trans = minus_log10, breaks = p_breaks, labels = as.character(p_breaks), sec.axis = sec_axis(~., breaks = q_limit, labels = paste("q =", fdr_limit)))
    } else {
      warning("None of the FDR-adjusted p-values are below the significance level, not plotting the horizontal line",
              call. = FALSE)
      pl <- pl +
        scale_y_continuous(trans = minus_log10, breaks = p_breaks,
                           labels = as.character(p_breaks))
    }

  } else {
    pl <- pl +
      scale_y_continuous(trans = minus_log10, breaks = p_breaks,
                         labels = as.character(p_breaks))
  }

  if (log2_x) {
    if (center_x_axis) {
      x_lim <- max(abs(log2(data[, effect_col])))
      x_lim <- c(2^(-x_lim), 2^x_lim)
    }
    pl <- pl +
      scale_x_continuous(trans = "log2", limits = x_lim)
  } else {
    if (center_x_axis) {
      x_lim <- max(abs(data[, effect_col]))
      x_lim <- c(-x_lim, x_lim)
    }
    pl <- pl +
      scale_x_continuous(limits = x_lim)
  }

  if (!is.null(label)) {
    if (label %in% colnames(data)) {
      label_data <- data[data[, label] < label_limit, ]
      pl <- pl +
        ggrepel::geom_label_repel(data = label_data,
                                  mapping = aes_string(label = label),
                                  seed = 313, alpha = 0.5, size = 3, force = 10,
                                  max.overlaps = 50
        )
    } else {
    warning("Label column not found, not plotting them")
    }
  }
  pl
}

# manhattan_plot(intestine_tse, x_col="Average_Mz", p_col="pvalues", q_col="qvalues", effect_dir = "sign_pvalues", color="pvalues", xlab="Average m/z", title=NULL)
setGeneric("manhattan_plot", signature = "object",
           function(object, x_col, p_col, effect_dir = NULL, q_col = NULL,
                    color = NULL, p_breaks = c(0.05, 0.01, 0.001, 1e-4), fdr_limit = 0.05, x_lim = NULL, y_lim = NULL,
                    color_scale = ggplot2::scale_color_viridis_c(),
                    title = "Manhattan plot", subtitle = NULL, ...) standardGeneric("manhattan_plot"))

setMethod("manhattan_plot", c(object = "SummarizedExperiment"),
          function(object, x_col, p_col, effect_dir = NULL, q_col = NULL,
                   color = NULL, p_breaks = c(0.05, 0.01, 0.001, 1e-4), fdr_limit = 0.05, x_lim = NULL, y_lim = NULL,
                   color_scale = ggplot2::scale_color_viridis_c(),
                   title = "Manhattan plot", subtitle = NULL, xlab=NULL, ...) {
            manhattan_plotter(as.data.frame(rowData(object)), x_col, p_col,
                              effect_dir, q_col, color, p_breaks, fdr_limit, x_lim, y_lim, color_scale, title, subtitle, xlab, ...)
          })

manhattan_plotter <- function(data, x_col, p_col, effect_dir, q_col, color,
                              p_breaks, fdr_limit, x_lim, y_lim, color_scale,
                              title, subtitle, xlab, ...) {

  if (min(data[, p_col]) > max(p_breaks)) {
    warning("All the p-values are larger than the p-value breaks supplied. Consider using larger p_breaks for plotting")
  }

  if (!is.null(effect_dir)) {
    data$y <- -log10(data[, p_col]) * sign(data[, effect_dir])
    p_labels <- outer(c(-1,1), p_breaks) %>% as.vector() %>% as.character()
    p_breaks <- outer(c(-1,1), -log10(p_breaks)) %>% as.vector()
    p_labels <- p_labels[order(p_breaks)]
    p_breaks <- sort(p_breaks)
  } else {
    data$y <- -log10(data[, p_col])
    p_labels <- as.character(p_breaks)
    p_breaks <- -log10(p_breaks)
    p_labels <- p_labels[order(p_breaks)]
    p_breaks <- sort(p_breaks)
  }

  pl <- ggplot(data, aes_string(x = x_col, y = "y", color = color)) +
    geom_point(...) +
    color_scale +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(),
          axis.ticks.y = element_blank()) +
    geom_hline(yintercept = 0, color = "grey") +
    labs(title = title, subtitle = subtitle, x=xlab, y = "Sign of effect * p-value")


  if(!is.null(q_col)) {

    if (any(data[, q_col] < fdr_limit)) {
      q_limit <- max(data[data[, q_col] < fdr_limit, p_col], na.rm = TRUE)
      if (!is.null(effect_dir)) {
        pl <- pl +
          geom_hline(yintercept = log10(q_limit), linetype = "dashed") +
          geom_hline(yintercept = -log10(q_limit), linetype = "dashed") +
          scale_y_continuous(breaks = p_breaks, labels = p_labels, limits = y_lim, sec.axis = sec_axis(~., breaks = c(log10(q_limit), -log10(q_limit)), labels = rep(paste("q =", fdr_limit), 2)))
      } else {
        pl <- pl +
          geom_hline(yintercept = -log10(q_limit), linetype = "dashed") +
          scale_y_continuous(breaks = p_breaks, labels = p_labels, limits = y_lim, sec.axis = sec_axis(~., breaks = -log10(q_limit), labels = paste("q =", fdr_limit)))
      }
    } else {
      warning("None of the FDR-adjusted p-values are below the significance level, not plotting the horizontal line",
              call. = FALSE)
      pl <- pl +
        scale_y_continuous(breaks = p_breaks, labels = p_labels, limits = y_lim)
    }
  } else {
    pl <- pl +
      scale_y_continuous(breaks = p_breaks, labels = p_labels, limits = y_lim)
  }
  pl
}

# cloud_plot(intestine_tse, p_col="pvalues", mz_col="Average_Mz", rt_col="Average_Rt_min_", color="pvalues")

setGeneric("cloud_plot", signature = "object",
           function(object, p_col = NULL, p_limit = NULL, mz_col = NULL,
                    rt_col = NULL, color = NULL, title = "m/z retention time", subtitle = NULL,
                    color_scale = ggplot2::scale_color_viridis_c(), ...) standardGeneric("cloud_plot"))

setMethod("cloud_plot", c(object = "SummarizedExperiment"),
          function(object, p_col = NULL, p_limit = NULL, mz_col = NULL,
                   rt_col = NULL, color = NULL,
                   title = "m/z vs retention time", subtitle = NULL,
                   color_scale = ggplot2::scale_color_viridis_c(),
                   all_features = FALSE) {
            cloud_plotter(as.data.frame(rowData(object)), p_col, p_limit, mz_col, rt_col, color, title, subtitle,
            color_scale, all_features)
          })

cloud_plotter <- function(x, p_col, p_limit, mz_col, rt_col, color, title,
                          subtitle, color_scale, all_features) {
  if (!is.null(p_limit) && !is.null(p_col)) {
    x <- x[x[, p_col] < p_limit, ]
    cat(paste("All features with p-values larger than", p_limit, "dropped.\n"))
  }
  if (is.null(mz_col) || is.null(rt_col)) {
    mz_rt_cols <- find_mz_rt_cols(x)
    mz_col <- mz_col %||% mz_rt_cols$mz_col
    rt_col <- rt_col %||% mz_rt_cols$rt_col
  }

  p <- ggplot(x, aes_string(x = rt_col, y = mz_col, size = p_col,
                            color = color)) +
    geom_point(alpha = 0.6) +
    scale_size_continuous(trans = minus_log10, range = c(0.5,3),
                          breaks = c(0.05, 0.01, 0.001, 1e-4),
                          labels = as.character(c(0.05, 0.01, 0.001, 1e-4))
    ) +
    theme_bw() +
    color_scale +
    labs(title = title, subtitle = subtitle,
         x = "Average RT (min)", y = "Average m/z", size = "p-value") +
    scale_x_continuous(breaks = seq(0, ceiling(max(x[, rt_col])))) +
    scale_y_continuous(breaks = seq(0, ceiling(max(x[, mz_col])), 50))

  if (length(unique(x$Split)) > 1) {
    cat("Multiple splits detected, plotting them to separate panes.\n")
    p <- p +
      facet_wrap(~Split, dir = "v")
  }
  p
}

# This function inputs a data.frame which must be constructed manually for interoperability with a variety of correlation functions for SE. For example with the phenomis package:
## intestine_tse <- hypotesting(intestine_tse, test.c="spearman", factor_names.vc="Injection_order", adjust_thresh.n = 1.00)

## injection_cor_df <- data.frame(x = rownames(intestine_tse), y = "Injection_order", Correlation_coefficient=rowData(intestine_tse)[, "spearman_Injection_order_cor"])

# effect_heatmap_plot(injection_cor_df, x_col = x, y = y, effect_col = "Correlation_coefficient")
effect_heatmap_plot <- function(data, x_col, y_col,
                                effect_col, p = NULL, p_limit = 0.1, point_size_range = c(1, 6),
                                log2_effect = FALSE, discretize_effect = FALSE, breaks = 5, reverse_y = TRUE,
                                use_coord_fixed = TRUE,
                                symmetric_aspect_ratio = TRUE, title = NULL, subtitle = NULL, fill_scale = NA) {

  if (!is.null(fill_scale)) {
    if (is.na(fill_scale)) {
      if (discretize_effect | class(data[, effect_col]) %in% c("factor", "character")) {
        fill_scale <- ggplot2::scale_fill_brewer(palette = "RdBu")
      } else {
        fill_scale <- ggplot2::scale_fill_distiller(palette = "RdBu")
      }
    }
  }

  legend_label <- effect_col
  if (log2_effect) {
    data[, effect_col] <- log2(data[, effect_col])
    legend_label <- paste0("log2(", effect_col, ")")
  }

  if (discretize_effect) {
    data[effect_col] <- cut(data[, effect_col], breaks = breaks, dig.lab = 2)
    data[effect_col] <- factor(data[, effect_col],
                               levels = rev(levels(data[, effect_col])))
  }

  ggp <- ggplot(data, aes_string(x = x_col , y = y_col, fill = effect_col)) +
    geom_tile() +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90)) +
    labs(x = "", y = "", fill = legend_label,
         title = title, subtitle = subtitle) +
    fill_scale

  if (use_coord_fixed) {
    ggp <- ggp +
      coord_fixed()
  }
  if (symmetric_aspect_ratio) {
    ggp <- ggp +
      theme(aspect.ratio=1)
  }

  if (!is.null(p)) {
    small_p <- data[data[, p] < p_limit, ]
    small_p[p] <- -log10(small_p[, p])
    ggp <- ggp +
      geom_point(aes_string(size = p), data = small_p,
                 colour = 'grey10', fill = 'grey30', alpha = 0.3, shape = 21) +
      scale_size(name = '-log10(p-value)', range = point_size_range)
  }
  ggp
}
