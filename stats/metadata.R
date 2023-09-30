###
# Analysis of the performance and metadata
#
# Outputs:
# 1. metadata.tex - sample demographics, average accuracy, longitudinal changes 
#    in accuracy, group differences in accuracy
# 2. fig1-dataset-overview.png - overview of the dataset
# 3. metadata.RData - intermediate results (accuracy ~ session LME, group 
#    differences t-test, mu power contrast t-values, data frames)
###


# Prefix for TeX output
prefix = "metadata"

# Load the data
data <- as.data.frame(read.csv(input_filename, na.strings = c('NaN', '')))
demo <- unique(data[,c('Subject', 'MBSRsubject', 'handedness', 'instrument', 
                       'athlete', 'handsport', 'hobby', 'gender', 'age')]) %>% 
  arrange(Subject)

# Load mu power contrasts
mu_contrast <- sanity$mu.power.diff
c(n_subjects, n_periods, n_sensors) %<-% dim(mu_contrast)

# Clear the file for exporting results to TeX
f <- file(output_filename, "w+")
close(f)


### Mu Power Contrasts
# Load one EEG dataset to import channel locations properly
EEG <- import_set(file.path(eeg.path, 'S1_Session_1.set'))
assert("Channel locations are complete", nrow(EEG$chan_info) == 60)

# Calculate t-values for each period and EEG channel
t.vals <- data.frame()
for (p in 1:n_periods) {
  for (ch in 1:n_sensors) {
    t.vals <- rbind(t.vals, data.frame(Period = p, 
      x = EEG$chan_info$x[ch], 
      y = EEG$chan_info$y[ch], 
      amplitude = t.test(mu_contrast[,p,ch])$statistic[['t']],
      label = EEG$chan_info$electrode[ch]))
  }
}
t.vals$Period = as.factor(t.vals$Period)
levels(t.vals$Period) <- c("rest", "prep", "fbend")


### Performance Metrics
# Look at the correlation between different columns
# All approaches for evaluating performance are more or less the same,
# proceed with percent valid correct then
corr_matrix <- cor(as.matrix(data[,1:16]))
names(data[,1:16])
p_corr <- corrplot(corr_matrix, tl.pos='n')

### Extract Accuracy Data
accuracy_data <- data[,c('Subject', 'Session', 'Accuracy', 'MBSRsubject')]
accuracy_data$Group <- factor(accuracy_data$MBSRsubject,
                              levels = c(0, 1),
                              labels = c('Control', 'MBSR'))
min_accuracy <- min(accuracy_data$Accuracy)

### Longitudinal Changes in Accuracy
fm <- scaleAndFitLME(accuracy_data, accuracy_vs_session_formula, 
                     columns_to_scale = c('Accuracy', 'Session'))
assert("Accuracy ~ Session did not converge",
       has_converged(fm))
coef <- lmerTest:::get_coefmat(fm)["Session",]
coef_ci <- confint(fm)
summary(fm)

### Difference in Performance between Groups
accuracy_avg <- accuracy_data %>%
  group_by(Subject) %>%
  reframe(Accuracy_avg = mean(Accuracy), 
          Subject = unique(Subject),
          Group = unique(Group))

mean_accuracy <- 100 * mean(accuracy_avg$Accuracy_avg)
mean_accuracy_mbsr <- 100 * mean(accuracy_avg$Accuracy_avg[accuracy_avg$Group == 'MBSR'])
mean_accuracy_control <- 100 * mean(accuracy_avg$Accuracy_avg[accuracy_avg$Group == 'Control'])

group_diff_ttest <- with(accuracy_avg, t.test(Accuracy_avg[Group == 'MBSR'],
                                              Accuracy_avg[Group == 'Control']))
group_diff_cohens_d <- cohens_d(accuracy_avg, Accuracy_avg ~ Group, 
                                comparisons = list(c("MBSR", "Control")),
                                var.equal = F)
print(group_diff_ttest)
print(group_diff_cohens_d)

### Mean Trial Length
triallength_avg <- data[,c('Subject', 'Session', 'MeanValidTrialLength')] %>%
  group_by(Subject) %>%
  summarise(TrialLength_avg = mean(MeanValidTrialLength))


### Group-average Accuracy for All Sessions
accuracy_session <- accuracy_data %>%
  group_by(Session) %>%
  summarise(Accuracy_avg = mean(Accuracy))
accuracy_first <- 100 * accuracy_session$Accuracy_avg[accuracy_session$Session == 1]
accuracy_last <- 100 * accuracy_session$Accuracy_avg[accuracy_session$Session == 11]


### Export the Results
tex.save(output_filename, "% Longitudinal Changes in Accuracy\n")
tex.save(output_filename, "AccuracySessionBeta", 
         format(as.numeric(coef[["Estimate"]]), digits = 2), prefix = prefix)
tex.save(output_filename, "AccuracySessionDf", 
         format(as.numeric(coef[["df"]]), digits = 3), prefix = prefix)
tex.save(output_filename, "AccuracySessionTvalue", 
         format(as.numeric(coef[["t value"]]), digits = 2), prefix = prefix)
tex.save(output_filename, "AccuracySessionPvalue", 
         format.pval(as.numeric(coef[["Pr(>|t|)"]]), 0.001, digits = 1), prefix = prefix)
tex.save(output_filename, "AccuracySessionBetaCIMin", 
         format(round(coef_ci["Session", 1], 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "AccuracySessionBetaCIMax", 
         format(round(coef_ci["Session", 2], 2), nsmall = 2), prefix = prefix)

tex.save(output_filename, "\n% Group Differences in Mean Accuracy\n")
tex.save(output_filename, "GroupDiffDf", 
         format(round(group_diff_ttest$parameter, 1), nsmall = 1), prefix = prefix)
tex.save(output_filename, "GroupDiffTvalue", 
         format(round(group_diff_ttest$statistic, 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "GroupDiffPvalue", 
         format.pval(round(group_diff_ttest$p.value, 2), 0.001, digits = 2), prefix = prefix)
tex.save(output_filename, "GroupDiffCIMin", 
         format(round(group_diff_ttest$conf.int[[1]], 2), digits = 2), prefix = prefix)
tex.save(output_filename, "GroupDiffCIMax", 
         format(round(group_diff_ttest$conf.int[[2]], 2), digits = 2), prefix = prefix)
# NOTE: order of groups in cohen's d and t-test does not match
tex.save(output_filename, "GroupDiffCohensD", 
         format(round(as.matrix(-group_diff_cohens_d[['effsize']])[[1,1]], 3), 
                digits = 1), prefix = prefix)

tex.save(output_filename, "\n% Mean Accuracy\n")
tex.save(output_filename, "MeanAccuracy", 
         format(round(mean_accuracy, 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "MeanAccuracyMBSR", 
         format(round(mean_accuracy_mbsr, 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "MeanAccuracyControl", 
         format(round(mean_accuracy_control, 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "MeanAccuracyFirst", 
         format(round(accuracy_first, 1), nsmall = 1), prefix = prefix)
tex.save(output_filename, "MeanAccuracyLast", 
         format(round(accuracy_last, 1), nsmall = 1), prefix = prefix)

tex.save(output_filename, "\n% Mean Trial Length\n")
tex.save(output_filename, "MeanTrialLength", 
         format(round(mean(triallength_avg$TrialLength_avg), 2), 
                nsmall = 2), 
         prefix = prefix)

### Plot the Results

# Histogram of Accuracy Across All Sessions
p_hist <- ggplot(data = accuracy_data, 
                 mapping = aes(x = Accuracy)) +
  geom_histogram(color = 'white') +
  ylab('Count') + 
  theme_classic()
  

# Longitudinal Changes in Average Accuracy
p_group <- accuracy_data %>%
  group_by(Session) %>%
  summarise(Accuracy_avg = mean(Accuracy), 
            Accuracy_se = sd(Accuracy) / sqrt(sum(!is.na(Accuracy))),
            Session = unique(Session)) %>%
  ungroup() %>%
  ggplot(data = .,
         aes(x = Session, y = Accuracy_avg)) +
  geom_pointrange(aes(ymin = Accuracy_avg - 1.96 * Accuracy_se, 
                      ymax = Accuracy_avg + 1.96 * Accuracy_se),
                  linewidth = 0.5, color = "black") +
  geom_smooth(method = "lm", se = F, color = 'blue', linewidth = 1) +
  scale_x_continuous(breaks = seq(1, 11, by = 2)) +
  scale_y_continuous(breaks = seq(0.6, 0.8, by = 0.1)) +
  ylab('Accuracy') +
  theme_classic()

# Load parts of the figure from other images
p_training <- ggdraw() +
  draw_image(file.path(asset.path, 'training-structure.png'),
             x = 0.05, width = 0.9)

p_trial <- ggdraw() +
  draw_image(file.path(asset.path, 'trial-structure.png'),
             x = 0.05, width = 0.9)

p_contrast <- ggplot(t.vals,
                     aes(x = x,
                         y = y,
                         fill = amplitude,
                         z = amplitude,
                         label = label)) +
  geom_topo(grid_res = 200,
            chan_size = rel(0.15), 
            head_size = rel(0.5),
            color = 'black',
            linetype = 'solid',
            linewidth = rel(0.1)) + 
  facet_wrap(. ~ Period, nrow = 1) +
  scale_fill_distiller(palette = "RdBu") + 
  theme_void() + 
  theme(legend.title = element_text(angle = -90), 
        legend.title.align = 0.5,
        strip.text.x = element_text(size = 10, margin = margin(b = 5))) +
  guides(fill = guide_colourbar(title.position = 'right',
                                title = 't-statistic')) +
  coord_equal()

# Raincloud plot - compare performance between groups
p_group_cmp <- ggplot(data = accuracy_avg,
                      aes(Group, Accuracy_avg, fill = Group)) +
  geom_rain(alpha = .5) +
  theme_classic() +
  scale_fill_brewer(palette = 'Dark2') +
  ylab('Accuracy') +
  guides(fill = 'none', color = 'none')

# Individual range of accuracy similar to Blankertz2010
shapes = c("Average" = 4, "Session" = 16)
labels = c("Average" = "Average Accuracy", "Session" = "Accuracy in Sessions")
p_individual <- accuracy_data %>%
  group_by(Subject) %>%
  reframe(Accuracy_avg = mean(Accuracy), 
          Accuracy_min = min(Accuracy),
          Accuracy_max = max(Accuracy),
          Accuracy = Accuracy,
          Subject = Subject) %>%
  ggplot(data = .,
         aes(x = reorder(Subject, Accuracy_avg), y = Accuracy_avg)) +
    geom_linerange(aes(ymin = Accuracy_min, ymax = Accuracy_max),
               linewidth = 0.25, color = "grey") + 
    geom_point(aes(y = Accuracy, shape = "Session"), size = 0.6, color = "grey") +
    geom_point(aes(shape = "Average"), size = 1) +
    xlab('Subject') + ylab('Accuracy') +
    expand_limits(y = c(min_accuracy, 1)) +
    scale_shape_manual(values = shapes, labels = labels) +
    theme_classic() +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      legend.position = c(0.8, 0.1),
      legend.title = element_blank()
    ) +
    guides(shape = guide_legend(override.aes = list(color = c('black', 'grey'))))



### Combine plots into dataset overview -- Figure 1
fig1_ABC <- plot_grid(p_training, p_trial, p_contrast,
                      ncol = 1, labels = c('A', 'B', 'C'),
                      align = 'h', axis = 'l')

fig1_DE <- plot_grid(p_group, p_group_cmp, nrow = 1, align = 'v',
                     labels = c('D', 'E'))
fig1_DEF <- plot_grid(fig1_DE, p_individual, ncol = 1, 
                      align = 'h', axis = 'l', 
                      rel_heights = c(0.33, 0.66), 
                      labels = c('', 'F'))
fig1 <- plot_grid(fig1_ABC, fig1_DEF, nrow = 1)
save_plot(file.path(plot.path, 'fig1-dataset-overview.pdf'), 
          fig1, base_width = 9, base_height = 6)
save_plot(file.path(plot.path, 'fig1-dataset-overview.png'), 
          fig1, bg = "white", base_width = 9, base_height = 6)


### Demographics

# Whole sample
numSubjects = nrow(demo)
numFemaleSubjects = nrow(demo[demo$gender == 'F',])
numRightHandedSubjects = nrow(demo[(demo$handedness == 'R') & !is.na(demo$handedness),])
meanAge = mean(demo$age, na.rm = T)
sdAge = sd(demo$age, na.rm = T)

# Export
tex.save(output_filename, "\n% Demographics - Whole Sample\n")
tex.save(output_filename, "numFemaleSubjects", 
         format(numFemaleSubjects), prefix = '')
tex.save(output_filename, "numRightHandedSubjects", 
         format(numRightHandedSubjects), prefix = '')
tex.save(output_filename, "meanAge", 
         format(round(meanAge, 1), nsmall = 1), prefix = '')
tex.save(output_filename, "sdAge", 
         format(round(sdAge, 1), nsmall = 1), prefix = '')


# MBSR group
demo_mbsr = demo[demo$MBSRsubject == 1,]
numMBSRSubjects = nrow(demo_mbsr)
numFemaleMBSRSubjects = nrow(demo_mbsr[demo_mbsr$gender == 'F',])
numRightHandedMBSRSubjects = nrow(demo_mbsr[(demo_mbsr$handedness == 'R') & !is.na(demo_mbsr$handedness),])
meanAgeMBSR = mean(demo_mbsr$age, na.rm = T)
sdAgeMBSR = sd(demo_mbsr$age, na.rm = T)

# Export
tex.save(output_filename, "\n% Demographics - MBSR Group\n")
tex.save(output_filename, "numMBSRSubjects", 
         format(numMBSRSubjects), prefix = '')
tex.save(output_filename, "numFemaleMBSRSubjects", 
         format(numFemaleMBSRSubjects), prefix = '')
tex.save(output_filename, "numRightHandedMBSRSubjects", 
         format(numRightHandedMBSRSubjects), prefix = '')
tex.save(output_filename, "meanAgeMBSR", 
         format(round(meanAgeMBSR, 1), nsmall = 1), prefix = '')
tex.save(output_filename, "sdAgeMBSR", 
         format(round(sdAgeMBSR, 1), nsmall = 1), prefix = '')

# Control group
demo_control = demo[demo$MBSRsubject == 0,]
numControlSubjects = nrow(demo_control)
numFemaleControlSubjects = nrow(demo_control[demo_control$gender == 'F',])
numRightHandedControlSubjects = nrow(demo_control[(demo_control$handedness == 'R') & !is.na(demo_control$handedness),])
meanAgeControl = mean(demo_control$age, na.rm = T)
sdAgeControl = sd(demo_control$age, na.rm = T)


# Export
tex.save(output_filename, "\n% Demographics - Control Group\n")
tex.save(output_filename, "numControlSubjects", 
         format(numControlSubjects), prefix = '')
tex.save(output_filename, "numFemaleControlSubjects", 
         format(numFemaleControlSubjects), prefix = '')
tex.save(output_filename, "numRightHandedControlSubjects", 
         format(numRightHandedControlSubjects), prefix = '')
tex.save(output_filename, "meanAgeControl", 
         format(round(meanAgeControl, 1), nsmall = 1), prefix = '')
tex.save(output_filename, "sdAgeControl", 
         format(round(sdAgeControl, 1), nsmall = 1), prefix = '')

# Save the results
save(t.vals, fm, coef_ci, group_diff_ttest, group_diff_cohens_d,
     accuracy_data, accuracy_avg, accuracy_session, demo, 
     file = file.path(r.path, 'metadata.RData'))