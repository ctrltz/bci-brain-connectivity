###
# Plotting functions
###


plotMultiverseBase <- function(df, val, sig, rows, cols, facet_rule,
                               val.name = 'Estimate', sig.name = 'Significant', 
                               lim = NULL, shapes = c("TRUE" = 16, "FALSE" = 1),
                               fontface = c('bold', 'plain'), strip.y.angle = 90) {
  # Plot the results of the multiverse analysis as a table for facilitating
  # visual comparison of results of different pipelines. Color codes the 
  # continuous value of interest (e.g., beta, t or correlation coefficient), 
  # dots indicate the significance of the effect of interest.
  #
  # The idea behind this function is to use nested rows and columns for coding
  # different pipeline steps. It is expected that the pipeline has at least
  # two steps (x and y), additional steps might be provided as a facet_rule
  #
  # Example: 
  #   pipeline: A -> B -> C -> D
  #   desired output: rows code B1, ..., Bn, grouped by A1, ..., An
  #                   columns code D1, ..., Dn, grouped by C1, ..., Cn
  #   arguments: x = D, y = B, facet_rule = A ~ C
  #
  # Parameters:
  #   df - data frame
  #   val - column to use for displaying the value of interest
  #   sig - column to use for displaying significance of the results
  #   rows, cols - pipeline steps (lowest level of hierarchy)
  #   facet_rule - pipeline steps (higher levels of hierarchy, 
  #                                cols3 + cols2 ~ rows3 + rows2)
  #   val.name - colorbar title for the value
  #   sig.name - legend title for the significance
  #   lim - data limit (calculated automatically if not provided)
  #   shapes - how to display significance
  #   fontface - optional styling for nested pipeline levels
  #   strip.y.angle - optional rotation of y labels
  # 
  # Returns:
  #   ggplot
  
  # Calculate data limits if not provided manually
  if (is.null(lim)) {
    lim = ceilingn(max(abs(df[,val])), 2)
  }
  
  # Calculate the number of terms in both sides of the facet rule
  facet_formula <- as.formula(facet_rule)
  parts <- unlist(strsplit(as.character(facet_formula), ' ~ '))
  n_lhs <- length(unlist(strsplit(parts[1], '[+]')))
  n_rhs <- length(unlist(strsplit(parts[2], '[+]')))
  
  # Prepare the main plot, colormap can be adapted separately
  p <- ggplot(df, aes(x = .data[[cols]], y = .data[[rows]])) +
    geom_raster(aes(fill = .data[[val]])) +
    geom_point(aes(shape = .data[[sig]]), color = "black",
               size = 3, fill = NA, stroke = 1) +
    scale_shape_manual(labels = c('TRUE' = 'Yes', 'FALSE' = 'No'),
                       limits = c('TRUE', 'FALSE'), values = shapes, drop = FALSE,
                       guide = guide_legend(ncol = 1, order = 2,
                                            title = sig.name, title.position = 'top')) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev) +
    guides(fill = guide_colorbar(order = 1, 
                                 title.position = 'top',
                                 title.hjust = 0.5)) +
    facet_nested(facet_formula, switch = "y",
                 strip = strip_nested(
                   text_x = elem_list_text(face = tail(fontface, n_rhs)),
                   text_y = elem_list_text(face = tail(fontface, n_lhs),
                                           angle = strip.y.angle),
                   by_layer_x = T, by_layer_y = T
                 )) +
    labs(fill = val.name) +
    theme_classic() +
    theme(
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom",
      panel.spacing = unit(0, "lines"),
      strip.background = element_blank(),
      strip.placement = "outside"
    )
  
  p
}


plotMultiverseSplit <- function(df, val, sig, facet_rule, 
                                val.name = bquote(beta / rho), lim = NULL) {
  # Plot the results of the split multiverse analysis
  # (setting up the appearance)
  #
  # Parameters:
  #   df - data frame
  #   val - column to use for displaying the value of interest
  #   sig - column to use for displaying significance of the results
  #   facet_rule - pipeline steps (higher levels of hierarchy, 
  #                                cols3 + cols2 ~ rows3 + rows2)
  #   val.name - legend title for the value of interest
  #   lim - data limit (calculated automatically if not provided)
  #
  # Returns:
  #   ggplot
  
  if (is.null(lim)) {
    lim = ceilingn(max(abs(df[,val])), 2)
  }
  
  p <- plotMultiverseBase(df, val, sig, 'Inverse', 'ROI_Method', 
                          facet_rule, val.name = val.name, lim = lim)
  
  p + scale_fill_gradient2(low = muted("blue"), high = muted("red"),
                           limits = c(-lim, lim), breaks = c(-lim, 0, lim)) +
    guides(y = guide_axis(angle = 90))
}


plotMultiverseJoint <- function(df, val, sig, facet_rule, lim = NULL) {
  # Plot the results of the joint multiverse analysis
  # (setting up the appearance)
  #
  # Parameters:
  #   df - data frame
  #   val - column to use for displaying the value of interest
  #   sig - column to use for displaying significance of the results
  #   facet_rule - pipeline steps (higher levels of hierarchy, 
  #                                cols3 + cols2 ~ rows3 + rows2)
  #   lim - data limit (calculated automatically if not provided)
  #
  # Returns:
  #   ggplot
  
  p <- plotMultiverseBase(df, val, sig, 'DummyY', 'DummyX', facet_rule, 
                          lim = lim, val.name = 'Consistency', strip.y.angle = 0)
  
  p + scale_fill_gradient2(low = "#FD8060", mid = "#FEE191", high = "#B0D8A4", 
                           midpoint = 0.5, limits = c(0, lim)) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
}


plotPipelineHistograms <- function(data, param, nrow, plot.path, prefix, width, height) {
  # Density plots of estimated values for different pipelines
  #
  # Parameters:
  #   data - data frame
  #   param - list with variable to plot (param$x) and filename (param$name)
  #   nrow - number of rows for facet wrap
  #   plot.path - where to save the plots
  #   prefix - prefix of the filename
  #   width, height - parameters for ggsave
  
  p <- ggplot(data, mapping = aes(x = .data[[param$x]])) +
    geom_histogram(mapping = aes(y = after_stat(density)),
                   bins = 50, alpha = 0.5) +
    geom_density() +
    facet_wrap(. ~ Pipeline, nrow = nrow) +
    plot_theme
  
  output_filename <- file.path(plot.path,
                               paste(prefix, param$name, sep = ''))
  ggsave(output_filename, p, width = width, height = height)
}


plotConnectivitySpectra <- function(df, measure, fdispmin = 3, fdispmax = 40, 
                                    linewidth = 0.7, x.axis = FALSE, y.axis = FALSE) {
  # Plot connectivity spectra for all pipelines and a specified measure
  #
  # Parameters:
  #   df - data frame with all connectivity results
  #   measure - which measure to plot (e.g., ImCoh_Within for within-hemisphere
  #     coherence)
  #   fdsipmin, fdispmax - limits (min/max, Hz) of frequencies to be displayed 
  #     in the spectra plots
  #   linewidth - line width to use
  #   x.axis, y.axis - whether to hide x/y axis labels (displayed by default)
  #
  # Returns:
  #   ggplot with spectra plotted for each of the pipelines
  
  colors <- c("eLORETA" = "#1984c5", "LCMV" = "#c23728")
  linetypes <- c("1SVD" = "solid", "3SVD" = "42", "AVG-F" = "11")
  
  freq_selection <- df$Freqs >= fdispmin & df$Freqs <= fdispmax
  
  p <- ggplot(df[freq_selection & df$Measure == measure,]) + 
    geom_line(aes(x = Freqs, y = Conn, color = Inverse, 
                  linetype = ROI_Method, group = Pipeline), 
              linewidth = linewidth) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = linetypes) +
    xlim(fdispmin, fdispmax) +
    xlab('Frequency (Hz)') + ylab('Phase Synchronization') +
    facet_nested(Type ~ Measure, scales = 'free', switch = 'y') +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.placement = 'outside',
          strip.text = element_text(face = 'bold'))
  
  if (!x.axis) {
    p <- p + 
      xlab('')
  }
  
  if (!y.axis) {
    p <- p + 
      theme(axis.title.y.left = element_blank(),
            strip.text.y = element_blank())
  }
  
  p
}


plotMultiverseScatterBetween <- function(df, x, y) {
  # Scatter plots for between-subject effects in all pipelines of the multiverse
  #
  # Parameters:
  #   df - data frame with results for all pipelines
  #   x, y - column names to use for x- and y-axis
  #
  # Returns:
  #   ggplot with multiple facets corresponding to all pipelines
  
  p_between <- ggplot(data = df) +
    geom_point(mapping = aes(x = .data[[x]], y = .data[[y]]), size = 0.5) +
    geom_line(mapping = aes(x = .data[[x]], y = .data[[y]]), color = "blue",
              stat="smooth", method = "lm", se = F) +
    facet_wrap(. ~ Pipeline, nrow = 4)
  p_between
}


plotMultiverseScatterWithin <- function(df, x, y) {
  # Scatter plots for within-subject effects in all pipelines of the multiverse
  #
  # Parameters:
  #   df - data frame with results for all pipelines
  #   x, y - column names to use for x- and y-axis
  #
  # Returns:
  #   ggplot with multiple facets corresponding to all pipelines
  
  p_within <- ggplot(data = df) +
    geom_point(mapping = aes(x = .data[[x]], y = .data[[y]], group = Subject, color = Subject), size = 0.5) +
    geom_line(mapping = aes(x = .data[[x]], y = .data[[y]], group = Subject), color = "black", alpha = 0.5, 
              stat="smooth", method = "lm", se = F, linewidth = 0.5) +
    geom_line(mapping = aes(x = .data[[x]], y = .data[[y]]), color = "blue",
              stat="smooth", method = "lm", se = F) +
    facet_wrap(. ~ Pipeline, nrow = 4) + 
    scale_colour_grey(start = 0.2, end = 0.8, aesthetics = "color", guide = "none")
  p_within
}


plotProcessingInteraction <- function(df, value_col, 
                                      group_col, g1, g2, 
                                      other_rule) {
  # Scatter plots for an interaction between two processing steps/methods
  # (comparing values in pairs of pipelines which differ only in one step)
  #
  # Parameters:
  #   df - data frame with results for all pipelines
  #   value_col - values from which column to use for plotting
  #   group_col - values from which column to use for splitting pipelines into 
  #     pairs
  #   g1, g2 - levels of group_col to contrast
  #   other_rule - how to arrange other processing steps in the facet grid (the
  #     logic is exactly the same as for facet_rule in plotMultiverseBase)
  #
  # Returns:
  #   ggplot with multiple facets corresponding to all pairs of pipelines
  
  g1_pipelines <- results[,group_col] == g1
  g2_pipelines <- results[,group_col] == g2
  
  cmp_df <- data.frame(
    Pipeline_g1 = df$Pipeline[g1_pipelines],
    Value_g1 = df[g1_pipelines, value_col],
    Pipeline_g2 = df$Pipeline[g2_pipelines],
    Value_g2 = df[g2_pipelines, value_col]
  )
  
  other_cols <- c('Inverse', 'Band', 'Mask', 'ROI_Method')
  for (other_col in other_cols) {
    if (other_col == group_col) {
      next
    }
    
    assert(paste("The values in column", other_col, "should be the same"),
           identical(df[g1_pipelines, other_col], 
                     df[g2_pipelines, other_col]))
    cmp_df[,other_col] <- df[g1_pipelines, other_col]
  }
  
  ggplot(cmp_df) +
    geom_abline(intercept = 0, slope = 1, color = 'blue', linetype = 'dashed') +
    geom_point(aes(x = Value_g1, y = Value_g2), color = 'darkgray', size = 0.5) +
    coord_fixed(ratio = 1) +
    xlab(paste(value_col, '-', g1)) + ylab(paste(value_col, '-', g2)) +
    facet_nested(other_rule, switch = 'y', 
                 labeller = labeller(Band = c('BB' = 'Broadband', 
                                              'NB' = 'Narrowband'))) + 
    theme_classic() + 
    theme(strip.background = element_blank(),
          strip.text = element_text(face = 'bold'),
          strip.placement = 'outside')
}


## Annotations based on p-values


mark <- function(pval, threshold = p.threshold) {
  # Choose bold font weight for significant p-values
  #
  # Parameters:
  #   pval - actual p-value
  #   threshold - significance threshold
  # 
  # Returns:
  #   "bold" if pval < threshold, "plain" otherwise
  
  if (pval < threshold) {
    "bold"
  } else {
    "plain"
  }
}


map_signif_level <- function(pval, signif_levels = c('*' = 0.05, 'n.s.' = 1)) {
  # Show significance with stars in plot annotations
  # 
  # Parameters:
  #   pval - actual pvalue
  #   signif_levels - named list with significance levels
  #
  # Returns:
  #   name of the first element in signif_levels among those greater than pval
  #
  # The implementation is copied from the source code of ggsignif
  
  names(which.min(signif_levels[which(signif_levels > pval)]))
}
