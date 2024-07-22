###
# Analysis of the performance and metadata
#
# Outputs:
# 1. metadata.tex - sample demographics, average accuracy, longitudinal changes 
#    in accuracy, group differences in accuracy, comparison of accuracy in 
#    different tasks
# 2. fig3-accuracy.png - overview of the participants' performance
# 3. figAsupp1-task-comparison.png - comparison of accuracy in different tasks
# 4. metadata.RData - intermediate results (accuracy ~ session LME, group 
#    differences t-test, mu power contrast t-values, data frames)
###


# Create output folder for images
prefix = "metadata"
output.folder <- file.path(plot.path, prefix)
dir.create(output.folder)

# Load the data
data <- as.data.frame(read.csv(input_filename, na.strings = c('NaN', '')))
demo <- unique(data[,c('Subject', 'MBSRsubject', 'handedness', 'instrument', 
                       'athlete', 'handsport', 'hobby', 'gender', 'age')]) %>% 
  arrange(Subject)

# Check for outliers in accuracy
assert("Checking for accuracy outliers", 
       sum(checkOutliers(data$Accuracy)) == 0)

# Clear the file for exporting results to TeX
f <- file(output_filename, "w+")
close(f)


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

beta0_acc <- lmerTest:::get_coefmat(fm)["(Intercept)", "Estimate"]
beta1_session <- lmerTest:::get_coefmat(fm)["Session", "Estimate"]

c(beta0_orig, beta1_orig) %<-% restoreBetas(beta0_acc, beta1_session,
                                            accuracy_data, "Session", "Accuracy")

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


### Class Balance 
class_balance <- data[,c('Subject', 'Session', 'ValidTrialsRight', 'ValidTrialsLeft')] %>%
  aggregate(. ~ Subject, mean)


### Group-average Accuracy for All Sessions
accuracy_session <- accuracy_data %>%
  group_by(Session) %>%
  summarise(Accuracy_avg = mean(Accuracy))
accuracy_first <- 100 * accuracy_session$Accuracy_avg[accuracy_session$Session == 1]
accuracy_last <- 100 * accuracy_session$Accuracy_avg[accuracy_session$Session == 11]


### Comparison of performance in different tasks
task_performance <- data[, c('Subject', 'Session', 'AccuracyTask1', 'AccuracyTask2', 'AccuracyTask3')] %>%
  rename(Task1 = AccuracyTask1, Task2 = AccuracyTask2, Task3 = AccuracyTask3)
task_avg_performance <- task_performance %>%
  aggregate(. ~ Subject, mean)


# Between-subject correlation of performance in tasks

task_cor_between_12 <- with(task_avg_performance, cor.test(Task1, Task2))
task_cor_between_13 <- with(task_avg_performance, cor.test(Task1, Task3))
task_cor_between_23 <- with(task_avg_performance, cor.test(Task2, Task3))

p_task_compare_12 <- ggplot(task_avg_performance,
                            aes(x = Task1, y = Task2)) +
  geom_point() +
  geom_smooth(method = 'lm', se = T, color = "blue") +
  annotate("text", x = 1, y = 0.425, 
           label = paste0('ρ = ', round(task_cor_between_12$estimate, 2)),
           hjust = 1, vjust = 0.5, size = 3.5,
           fontface = mark(task_cor_between_12$p.value)) +
  xlab('Accuracy (task 1)') + ylab('Accuracy (task 2)') +
  ggtitle('Task 1 ~ Task 2') +
  expand_limits(x = c(0.4, 1), y = c(0.4, 1)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 10))


p_task_compare_13 <- ggplot(task_avg_performance,
                            aes(x = Task1, y = Task3)) +
  geom_point() +
  geom_smooth(method = 'lm', se = T, color = "blue") +
  annotate("text", x = 1, y = 0.225, 
           label = paste0('ρ = ', round(task_cor_between_13$estimate, 2)),
           hjust = 1, vjust = 0.5, size = 3.5,
           fontface = mark(task_cor_between_13$p.value)) +
  xlab('Accuracy (task 1)') + ylab('Accuracy (task 3)') +
  ggtitle('Task 1 ~ Task 3') +
  expand_limits(x = c(0.4, 1), y = c(0.2, 1)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 10))


p_task_compare_23 <- ggplot(task_avg_performance,
                            aes(x = Task2, y = Task3)) +
  geom_point() +
  geom_smooth(method = 'lm', se = T, color = "blue") +
  annotate("text", x = 1, y = 0.225, 
           label = paste0('ρ = ', round(task_cor_between_23$estimate, 2)),
           hjust = 1, vjust = 0.5, size = 3.5,
           fontface = mark(task_cor_between_23$p.value)) +
  xlab('Accuracy (task 2)') + ylab('Accuracy (task 3)') +
  ggtitle('Task 2 ~ Task 3') +
  expand_limits(x = c(0.4, 1), y = c(0.2, 1)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 10))


## Within-subject correlation of performance in tasks

task_corr <- task_performance %>%
  nest_by(Subject) %>%
  mutate(corr_12 = with(data, cor(Task1, Task2)),
         corr_13 = with(data, cor(Task1, Task3)),
         corr_23 = with(data, cor(Task2, Task3)))
n_subjects <- nrow(task_corr)

within_corr_df <- data.frame(
  Tasks = c(
    rep('Task 1 ~ Task 2', n_subjects),
    rep('Task 1 ~ Task 3', n_subjects),
    rep('Task 2 ~ Task 3', n_subjects)
  ),
  Subject = rep(task_corr$Subject, 3),
  Correlation = unlist(c(task_corr$corr_12, task_corr$corr_13, task_corr$corr_23))
)

median_within_corr_df <- within_corr_df %>% 
  aggregate(Correlation ~ Tasks, median)
c(within_corr_12, within_corr_13, within_corr_23) %<-% median_within_corr_df$Correlation

p_task_compare_within <- ggplot(within_corr_df, aes(x = 1, y = Correlation, fill = Tasks)) +
  geom_rain(alpha = .5) +
  facet_wrap(. ~ Tasks, nrow = 1) +
  xlab('') + ylab('Within-subject correlation') +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x = element_text(size = 10)) +
  guides(fill = 'none', color = 'none') + 
  coord_flip()



fig_task_cmp_between <- plot_grid(p_task_compare_12, p_task_compare_13, p_task_compare_23,
                                  nrow = 1, align = 'v')
fig_task_cmp <- plot_grid(fig_task_cmp_between, p_task_compare_within, ncol = 1,
                          align = 'h', axis = 'l', labels = c('A', 'B'))

save_plot(file.path(output.folder, 'figAsupp1-task-comparison.png'),
          fig_task_cmp, base_width = 9, base_height = 6, bg = "white")


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
# NOTE: order of groups in cohen's d and t-test does not match, multiplying the
# t-value by minus one
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

tex.save(output_filename, "\n% Class Balance\n")
tex.save(output_filename, "MeanTrialsRight", 
         format(round(mean(class_balance$ValidTrialsRight), 2), 
                nsmall = 2), 
         prefix = prefix)
tex.save(output_filename, "SDTrialsRight", 
         format(round(sd(class_balance$ValidTrialsRight), 2), 
                digits = 2), 
         prefix = prefix)
tex.save(output_filename, "MeanTrialsLeft", 
         format(round(mean(class_balance$ValidTrialsLeft), 2), 
                nsmall = 2), 
         prefix = prefix)
tex.save(output_filename, "SDTrialsLeft", 
         format(round(sd(class_balance$ValidTrialsLeft), 2), 
                nsmall = 2), 
         prefix = prefix)

tex.save(output_filename, "\n% Performance in Different Tasks\n")
tex.save(output_filename, "CorrFirstSecondTask", 
         format(round(task_cor_between_12$estimate, 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "CorrFirstSecondTaskCIMin", 
         format(round(task_cor_between_12$conf.int[[1]], 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "CorrFirstSecondTaskCIMax", 
         format(round(task_cor_between_12$conf.int[[2]], 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "CorrFirstThirdTask", 
         format(round(task_cor_between_13$estimate, 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "CorrFirstThirdTaskCIMin", 
         format(round(task_cor_between_13$conf.int[[1]], 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "CorrFirstThirdTaskCIMax", 
         format(round(task_cor_between_13$conf.int[[2]], 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "CorrSecondThirdTask", 
         format(round(task_cor_between_23$estimate, 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "CorrSecondThirdTaskCIMin", 
         format(round(task_cor_between_23$conf.int[[1]], 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "CorrSecondThirdTaskCIMax", 
         format(round(task_cor_between_23$conf.int[[2]], 2), nsmall = 2), prefix = prefix)

tex.save(output_filename, "CorrWithinFirstSecondTask", 
         format(round(within_corr_12, 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "CorrWithinFirstThirdTask", 
         format(round(within_corr_13, 2), nsmall = 2), prefix = prefix)
tex.save(output_filename, "CorrWithinSecondThirdTask", 
         format(round(within_corr_23, 2), nsmall = 2), prefix = prefix)

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
  expand_limits(y = c(0.6, 0.8)) +
  ylab('Accuracy') +
  theme_classic()


# Raincloud plot - compare performance between groups
p_group_cmp <- ggplot(data = accuracy_avg,
                      aes(Group, Accuracy_avg, fill = Group)) +
  geom_rain(alpha = .5) +
  geom_signif(comparisons = list(c("Control", "MBSR")), y_position = 0.95,
              tip_length = 0, test = "t.test", textsize = 4,
              map_signif_level = c('*' = 0.05, 'n.s.' = 1)) +
  expand_limits(y = c(0.5, 1)) +
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
    expand_limits(y = c(0.2, 1)) +
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
fig3_AB <- plot_grid(p_group, p_group_cmp, ncol = 1, 
                     align = 'h', axis = 'l',
                     labels = c('A', 'B'))
fig3 <- plot_grid(fig3_AB, p_individual, 
                  nrow = 1, rel_widths = c(1, 2),
                  labels = c('', 'C'))
save_plot(file.path(output.folder, 'fig3-accuracy.pdf'), 
          fig3, base_width = 9, base_height = 4)
save_plot(file.path(output.folder, 'fig3-accuracy.png'), 
          fig3, bg = "white", base_width = 7, base_height = 4)

fig3_parts = list(
  group = p_group, group_cmp = p_group_cmp, individual = p_individual
)


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
save(fm, coef_ci, group_diff_ttest, group_diff_cohens_d,
     accuracy_data, accuracy_avg, accuracy_session, demo, 
     task_cor_between_12, task_cor_between_13, task_cor_between_23,
     within_corr_df, file = file.path(r.path, 'metadata.RData'))