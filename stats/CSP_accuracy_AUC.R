###
# Comparison of online accuracy and offline (CSP-based) accuracy and AUC
#
# Outputs:
# 1. CSP_accuracy_AUC.tex - comparison of offline performance metrics with 
#    chance levels in all time windows, correlations with online accuracy
#    (between- and within-subject)
# 2. figBsupp1-csp-accuracy-auc.png - comparison of accuracy in different tasks
# 3. CSP_accuracy_AUC.RData - intermediate results (group t-test table,
#    correlations)
###

# Load the data
data <- readMat(input_filename)
csp_df <- as.data.frame(data$csp.accuracy.data)
colnames(csp_df) <- unlist(data$labels)
csp_df$Period <- as.factor(csp_df$Period)
period_names <- c('rest', 'target', 'feedback')
levels(csp_df$Period) <- period_names

# Create output folder for images
prefix = "csp"
output.folder <- file.path(plot.path, prefix)
dir.create(output.folder)

# Clear the file for exporting results to TeX
f <- file(output_filename, "w+")
close(f)


csp_df_avg <- csp_df %>%
  aggregate(. ~ Subject + Period, mean)


### Statistics

# Accuracy and AUC in different time windows
n.csp.mc <- 3           # 3 time windows considered
chance.level <- 0.5     # for both accuracy and AUC
csp_group_ttest <- csp_df_avg %>%
  nest_by(Period) %>%
  mutate(test.acc = list(with(data, t.test(Accuracy, rep(0.5, 62)))),
         test.auc = list(with(data, t.test(AUC, rep(0.5, 62))))) %>%
  mutate(mean.acc = test.acc$estimate[[1]],
         t.acc = test.acc$statistic,
         df.acc = test.acc$parameter,
         p.acc = test.acc$p.value,
         ci.min.acc = test.acc$conf.int[[1]] + chance.level,
         ci.max.acc = test.acc$conf.int[[2]] + chance.level,
         mean.auc = test.auc$estimate[[1]],
         t.auc = test.auc$statistic,
         df.auc = test.auc$parameter,
         p.auc = test.auc$p.value,
         ci.min.auc = test.auc$conf.int[[1]] + chance.level,
         ci.max.auc = test.auc$conf.int[[2]] + chance.level) %>%
  select(!c(data, test.acc, test.auc))


# Use the feedback interval for comparing online and offline accuracy

csp_df_fbend <- csp_df_avg[csp_df_avg$Period == period_names[[3]],]

csp_corr_within <- csp_df[csp_df$Period == period_names[[3]],] %>% 
  nest_by(Subject) %>%
  mutate(cor.acc = list(with(data, cor.test(Accuracy, Online))),
         cor.auc = list(with(data, cor.test(AUC, Online)))) %>%
  mutate(c.acc = cor.acc$estimate, 
         c.auc = cor.auc$estimate)


# Between-subject correlation on average values of accuracy/AUC

corr_csp_acc_between <- with(csp_df_fbend, cor.test(Accuracy, Online))
corr_csp_acc_between

corr_csp_auc_between <- with(csp_df_fbend, cor.test(AUC, Online))
corr_csp_auc_between


### Prepare the output table
csp_group_ttest_acc <- csp_group_ttest %>%
  select(Period, mean.acc, t.acc, df.acc, p.acc, ci.min.acc, ci.max.acc) %>%
  mutate(p.mc = p.acc * n.csp.mc) %>%
  rename(Mean = mean.acc, "t-value" = t.acc,
         df = df.acc, "p-value" = p.acc,
         "p-value (corr.)" = p.mc,
         "CIMin" = ci.min.acc, "CIMax" = ci.max.acc)

csp_group_ttest_auc <- csp_group_ttest %>%
  select(Period, mean.auc, t.auc, df.auc, p.auc, ci.min.auc, ci.max.auc) %>%
  mutate(p.mc = p.auc * n.csp.mc) %>%
  rename(Mean = mean.auc, "t-value" = t.auc,
         df = df.auc, "p-value" = p.auc,
         "p-value (corr.)" = p.mc,
         "CIMin" = ci.min.auc, "CIMax" = ci.max.auc)

csp_group_ttest_table <- bind_rows(csp_group_ttest_acc, csp_group_ttest_auc) %>%
  mutate_at(vars("t-value", "Mean", "CIMin", "CIMax"), round, 2) %>%
  mutate_at(vars("p-value", "p-value (corr.)"), round, 3) %>%
  mutate(Significant = `p-value` < 0.05,
         Significant.MC = `p-value (corr.)` < 0.05,
         star = if_else(Significant.MC, '*', ''))

# Add the metric column
csp_group_ttest_table$Metric = c(
  rep('Accuracy', 3), rep('AUC', 3)
)

# Format the confidence intervals
csp_group_ttest_table$CIDisp <- paste(
  "$[", format(csp_group_ttest_table$CIMin, nsmall = 2), ",", 
  format(csp_group_ttest_table$CIMax, nsmall = 2), "]$", sep = "")

# Format the p-values: bold denotes significance BEFORE correction for multiple
# comparison, stars - AFTER the correction
csp_group_ttest_table$`p-value` <- with(csp_group_ttest_table,
                                        cell_spec(paste(
                                          if_else(`p-value` < 0.001, 
                                                  '$<$0.001', 
                                                  format(`p-value`, digits = 3)), 
                                          star, 
                                          sep = ''), 
                                          format = 'latex', 
                                          bold = Significant, 
                                          escape = F))


### Export the Results
tex.save(output_filename, "% Group-level t-tests for accuracy and AUC")
csp_group_ttest_latex <- csp_group_ttest_table %>%
  rename("95\\% CI" = "CIDisp") %>%
  select(Metric, Period, Mean, "t-value", df, "p-value", "95\\% CI") %>%
  knitr::kable(format = "latex", booktabs = T, escape = F, align = 'c',
               linesep = c("", "", "", "", "", "\\midrule")) %>%
  collapse_rows(columns = 1) %>%
  gsub("\\cmidrule{2-7}\n", "", ., fixed = T) %>%
  tex.save(output_filename, "GroupTTestSummary", ., prefix = prefix)

tex.save(output_filename, "\n% Between-Subject Correlation\n")
tex.save(output_filename, "CorrAccuracyBetween", 
         format(as.numeric(corr_csp_acc_between$estimate), digits = 2), prefix = prefix)
tex.save(output_filename, "CorrAccuracyBetweenCIMin", 
         format(as.numeric(corr_csp_acc_between$conf.int[[1]]), digits = 2), prefix = prefix)
tex.save(output_filename, "CorrAccuracyBetweenCIMax", 
         format(as.numeric(corr_csp_acc_between$conf.int[[2]]), digits = 2), prefix = prefix)
tex.save(output_filename, "CorrAUCBetween", 
         format(as.numeric(corr_csp_auc_between$estimate), digits = 2), prefix = prefix)
tex.save(output_filename, "CorrAUCBetweenCIMin", 
         format(as.numeric(corr_csp_auc_between$conf.int[[1]]), digits = 2), prefix = prefix)
tex.save(output_filename, "CorrAUCBetweenCIMax", 
         format(as.numeric(corr_csp_auc_between$conf.int[[2]]), digits = 2), prefix = prefix)


tex.save(output_filename, "\n% Within-Subject Correlation\n")
tex.save(output_filename, "CorrAccuracyWithin", 
         format(as.numeric(median(csp_corr_within$c.acc)), digits = 2), prefix = prefix)
tex.save(output_filename, "CorrAUCWithin", 
         format(as.numeric(median(csp_corr_within$c.auc)), digits = 2), prefix = prefix)


### Plots

# Average accuracy and AUC for each time period

acc_signif_annot <- lapply(period_names,
                           \(x) map_signif_level(csp_group_ttest_acc$`p-value (corr.)`[csp_group_ttest_acc$Period == x]))
p_acc_periods <- ggplot(csp_df_avg, aes(x = Period, y = Accuracy, fill = Period)) + 
  geom_abline(intercept = 0.5, slope = 0,
              color = "#aaaaaa", linetype = "dashed") +
  geom_rain(alpha = .5) + 
  annotate("text", x = period_names, y = 1, 
           label = acc_signif_annot,
           vjust = 0.25, size = 4) +
  expand_limits(y = c(0.4, 1)) +
  xlab('Time window') + ylab('Accuracy (offline)') +
  theme_classic() +
  guides(fill = 'none')

auc_signif_annot <- lapply(period_names,
                           \(x) map_signif_level(csp_group_ttest_auc$`p-value (corr.)`[csp_group_ttest_auc$Period == x]))
p_auc_periods <- ggplot(csp_df_avg, aes(x = Period, y = AUC, fill = Period)) + 
  geom_abline(intercept = 0.5, slope = 0,
              color = "#aaaaaa", linetype = "dashed") +
  geom_rain(alpha = .5,
            point.args = list(shape = 24, fill = "black", alpha = .5)) +
  annotate("text", x = period_names, y = 1, 
           label = auc_signif_annot,
           vjust = 0.25, size = 4) +
  expand_limits(y = c(0.4, 1)) +
  xlab('Time window') + ylab('AUC (offline)') +
  theme_classic() + 
  guides(fill = 'none')


# Between-subject correlation of online and offline accuracy/AUC

fbend_color <- hue_pal()(3)[[3]]

p_acc_between <- ggplot(csp_df_fbend, aes(x = Online, y = Accuracy)) +
  geom_abline(intercept = 0, slope = 1,
              color = "#aaaaaa", linetype = "dashed") +
  geom_point(color = fbend_color) + 
  annotate("text", x = 1, y = 0.425, 
           label = paste0('ρ = ', round(corr_csp_acc_between$estimate, 2)),
           hjust = 1, vjust = 0.5, size = 3.5,
           fontface = mark(corr_csp_acc_between$p.value)) +
  expand_limits(x = c(0.4, 1), y = c(0.4, 1)) +
  xlab('Accuracy (online)') + ylab('Accuracy (offline)') +
  theme_classic()

p_auc_between <- ggplot(csp_df_fbend, aes(x = Online, y = AUC)) +
  geom_abline(intercept = 0, slope = 1,
              color = "#aaaaaa", linetype = "dashed") +
  geom_point(color = fbend_color, fill = fbend_color, shape = 24) + 
  annotate("text", x = 1, y = 0.425, 
           label = paste0('ρ = ', round(corr_csp_auc_between$estimate, 2)),
           hjust = 1, vjust = 0.5, size = 3.5,
           fontface = mark(corr_csp_auc_between$p.value)) +
  expand_limits(x = c(0.4, 1), y = c(0.4, 1)) +
  xlab('Accuracy (online)') + ylab('AUC (offline)') +
  theme_classic()


# Within-subject ICC of online and offline accuracy/AUC

p_acc_within <- ggplot(csp_corr_within, aes(x = 1, y = c.acc)) + 
  geom_rain(alpha = .5, fill = fbend_color) +
  ylab('Within-subject correlation of\n offline and online accuracy') +
  coord_flip() + 
  theme_classic() + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_auc_within <- ggplot(csp_corr_within, aes(x = 1, y = c.auc)) + 
  geom_rain(alpha = .5, fill = fbend_color,
            point.args = list(shape = 24, fill = "black", alpha = .5)) +
  ylab('Within-subject correlation of\n AUC and online accuracy') +
  coord_flip() + 
  theme_classic() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank())


# Combine the plots
p_label_acc = ggdraw() +
  draw_label("Offline accuracy (CSP + rLDA)",
             fontface = 'bold', y = 0.75)
p_label_auc = ggdraw() +
  draw_label("Offline AUC (CSP + rLDA)",
             fontface = 'bold', y = 0.25, vjust = 0)


fig_csp_accuracy <- plot_grid(p_acc_periods, p_acc_between, p_acc_within,
                              nrow = 1, align = 'hv', labels = c('A', 'B', 'C'))

fig_csp_auc <- plot_grid(p_auc_periods, p_auc_between, p_auc_within, 
                         nrow = 1, align = 'hv', labels = c('D', 'E', 'F'))

fig_csp <- plot_grid(p_label_acc, fig_csp_accuracy, p_label_auc, fig_csp_auc,
                     ncol = 1, rel_heights = c(0.2, 1, 0.2, 1), 
                     align = 'h', axis = 'l')

save_plot(file.path(output.folder, 'figBsupp1-csp-accuracy-auc.png'),
          fig_csp, base_width = 9, base_height = 7, bg = "white")


# Save the results
save(csp_group_ttest, csp_corr_within, corr_csp_acc_between,
     file = file.path(r.path, 'CSP_accuracy_AUC.RData'))