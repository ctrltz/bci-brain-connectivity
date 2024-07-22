%%%
% Sanity checks related to power and connectivity, data for plots in R
%
% Outputs:
% 1. BCI_MI_sanity_checks.mat with the following variables:
%
%   demo_snr - PSD and FOOOF fit for SNR panel in the pipeline overview
%   example_snr - Laplace PSD and SNR for three exemplary subjects
%   mu_power_diff - mu power contrast between imaginary movements
%   conn_data, conn_labels, pipeline_labels - within- and across-hemisphere
%       connectivity spectra
%%%

%% Mu Power Contrasts - Load the Data
load([savedata 'BCI_MI_mu_power_contrasts.mat']);

%% Mu Power Contrasts - Sensor Space
n_analyzed_sessions = sum(results_mu.analyzed_sessions, [1 2]);

% take analyzed and good sessions into account
spec_sensors = results_mu.spec_sensors;
for subject = 1:n_subjects
    for session = 1:n_sessions
        to_keep = results_mu.analyzed_sessions(subject, session) && ...
            good_sessions(subject, session);

        if ~to_keep
            disp([subject session]);
            spec_sensors(subject, session, :, :, :, :) = NaN;
        end
    end
end
spec_sensors_avg = squeeze(mean(spec_sensors, 2, 'omitnan'));
assert(all(~isnan(spec_sensors_avg), 'all'));

mu_power = squeeze(mean(spec_sensors_avg(:, :, :, :, mu_band_bins), 5));
mu_power_diff = squeeze(mu_power(:, :, 1, :) - mu_power(:, :, 2, :));

tvals = zeros(n_periods, n_sensors);
pvals = zeros(n_periods, n_sensors);
for p = 1:n_periods
    [~, pvals(p, :), ~, stats] = ttest(squeeze(mu_power(:, p, 1, :) - mu_power(:, p, 2, :)));
    tvals(p, :) = stats.tstat;
end
lim = max(abs(tvals), [], 'all');

for p = 1:n_periods
    h = figure('Position', [840 666 605 220]);
    topoplot(tvals(p, :), all_chanlocs, 'maplimits', [-lim lim]);
    ax = colorbar;
    ax.Label.String = 't-statistic';
    ax.Label.Position(1) = 3; 
    ax.Label.FontSize = 10;
    ax.Label.Rotation = 270;
end

%% Mu Power Contrasts - Clear the Data
clear results_mu;

%% Laplace - Load the Data
load([savedata 'BCI_MI_Laplace_SNR.mat']);

%% PSD Example with fooof fit for pipeline overview
spec_demo = squeeze(mean(results_lap.spec_laplace(:, 1, :, :), [1 2 3]));
fooof_results = fooof(freqs, spec_demo, [3 40], fooof_settings, return_model);
spec_orig = 10 .^ fooof_results.power_spectrum;
spec_noise = 10 .^ fooof_results.ap_fit;
fit_freqs = fooof_results.freqs;
mu_bins = fit_freqs >= mu_band(1) & fit_freqs <= mu_band(2);
mu_freqs = fit_freqs(mu_bins);

figure; hold on;
p1 = plot(fit_freqs, spec_orig, 'LineWidth', 2, 'Color', 'k');
p2 = plot(fit_freqs, spec_noise, 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
fill([mu_freqs fliplr(mu_freqs)], [spec_noise(mu_bins) fliplr(spec_orig(mu_bins))], 'b');
fill([mu_freqs fliplr(mu_freqs)], [zeros(size(mu_freqs)) fliplr(spec_noise(mu_bins))], 'r');
legend([p1 p2], {'PSD', '1/f fit'}, 'Location', 'northeast', 'Box', 'off');
xlim([3 40]);
xlabel('Frequency (Hz)'); ylabel('PSD (\muV^2/Hz)');

demo_snr = struct();
demo_snr.spec_orig = spec_orig;
demo_snr.spec_noise = spec_noise;
demo_snr.fit_freqs = fit_freqs;
demo_snr.mu_freqs = mu_freqs;
demo_snr.mu_bins = mu_bins;

%% Laplace SNR (C3/C4-Laplace, resting-state)
spec_laplace_avg = squeeze(mean(results_lap.spec_laplace, [2 3]));
snr_laplace_avg = squeeze(mean(10*log10(results_lap.snr_laplace), ...
    [2 3], 'omitnan'));
figure; hold on;
plot(freqs, 10*log10(spec_laplace_avg(subj_low_snr, :)), ...
    'LineWidth', 2, 'DisplayName', sprintf('SNR = %.2f', snr_laplace_avg(subj_low_snr)));
plot(freqs, 10*log10(spec_laplace_avg(subj_med_snr, :)), ...
    'LineWidth', 2, 'DisplayName', sprintf('SNR = %.2f', snr_laplace_avg(subj_med_snr)));
plot(freqs, 10*log10(spec_laplace_avg(subj_high_snr, :)), ...
    'LineWidth', 2, 'DisplayName', sprintf('SNR = %.2f', snr_laplace_avg(subj_high_snr)));
xlabel('Frequency (Hz)'); ylabel('10*log10(PSD)');
xlim([1 40]);
legend('Location', 'northeast');

example_snr = struct();
example_snr.freqs = freqs;
example_snr.spec_avg = spec_laplace_avg;
example_snr.snr_avg = snr_laplace_avg;
example_snr.subj_low_snr = subj_low_snr;
example_snr.subj_med_snr = subj_med_snr;
example_snr.subj_high_snr = subj_high_snr;

%% Laplace - Clear to Save RAM
clear results_lap;

%% Multiverse - Load the Data
load([savedata 'BCI_MI_spectra_rois_multiverse.mat']);
load([savedata 'BCI_MI_connectivity_multiverse.mat']);

%% Multiverse analysis
analyzed_sessions = results_conn.analyzed_sessions(:) & good_sessions(:);
bb_pipelines = cellfun(@(x) contains(x, 'bb'), multiverse_labels);
bb_pipeline_idx = find(bb_pipelines);
n_bb_pipelines = sum(bb_pipelines);

%% Power spectra
spec_rois_packed = reshape(cellfun(@(x) mean(x, 1)', ...
    squeeze(results_spec.spec_rois(bb_pipelines, :, :)), 'UniformOutput', false),...
    n_bb_pipelines, []);
spec_rois_packed = spec_rois_packed(:, analyzed_sessions > 0);
spec_rois = zeros([size(spec_rois_packed) n_freqs]);
for i = 1:n_bb_pipelines
    spec_rois(i, :, :) = cell2mat(squeeze(spec_rois_packed(i, :)))';
end
spec_rois = squeeze(mean(spec_rois, 2));

figure;
for i = 1:n_bb_pipelines
    subplot(4, 3, i); hold on;
    p1 = plot(freqs, 10*log10(spec_rois(i, :)), 'k');
    title(multiverse_labels{bb_pipeline_idx(i)});
end
sgtitle('10*log10(PSD)');

%% Connectivity spectra
% abs(ImCoh)
absicoh_freq_rest = results_conn.absicoh_freq(bb_pipelines, :, :, :, :, :);
within_absicoh_freq_rest = reshape((absicoh_freq_rest(:, :, :, :, 1, 2) + ...
    absicoh_freq_rest(:, :, :, :, 3, 4)) / 2, n_bb_pipelines, [], n_freqs);
within_absicoh_freq_rest = squeeze(mean(within_absicoh_freq_rest(:, analyzed_sessions > 0, :), 2))';
across_absicoh_freq_rest = reshape((absicoh_freq_rest(:, :, :, :, 1, 3) + ...
    absicoh_freq_rest(:, :, :, :, 1, 4) + absicoh_freq_rest(:, :, :, :, 2, 3) + ...
    absicoh_freq_rest(:, :, :, :, 2, 4)) / 4, n_bb_pipelines, [], n_freqs);
across_absicoh_freq_rest = squeeze(mean(across_absicoh_freq_rest(:, analyzed_sessions > 0, :), 2))';

figure;
for i = 1:n_bb_pipelines
    subplot(5, 3, i); hold on;
    p1 = plot(freqs, within_absicoh_freq_rest(:, i), 'b');
    p2 = plot(freqs, across_absicoh_freq_rest(:, i), 'r');
    title(multiverse_labels{bb_pipeline_idx(i)});
end
ax = subplot(5, 3, 14);
ax.Visible = 'off';
lgd = legend(ax, [p1 p2], {'Within-hemisphere', 'Across-hemisphere'}, ...
    'Location', 'northoutside');
title(lgd, 'Connections');
sgtitle('abs(ImCoh)');

% abs(Coh)
coh_freq_rest = results_conn.coh_freq(bb_pipelines, :, :, :, :, :);
within_coh_freq_rest = reshape((coh_freq_rest(:, :, :, :, 1, 2) + ...
    coh_freq_rest(:, :, :, :, 3, 4)) / 2, n_bb_pipelines, [], n_freqs);
within_coh_freq_rest = squeeze(mean(within_coh_freq_rest(:, analyzed_sessions > 0, :), 2))';
across_coh_freq_rest = reshape((coh_freq_rest(:, :, :, :, 1, 3) + ...
    coh_freq_rest(:, :, :, :, 1, 4) + coh_freq_rest(:, :, :, :, 2, 3) + ...
    coh_freq_rest(:, :, :, :, 2, 4)) / 4, n_bb_pipelines, [], n_freqs);
across_coh_freq_rest = squeeze(mean(across_coh_freq_rest(:, analyzed_sessions > 0, :), 2))';

figure;
for i = 1:n_bb_pipelines
    subplot(5, 3, i); hold on;
    p1 = plot(freqs, within_coh_freq_rest(:, i), 'b');
    p2 = plot(freqs, across_coh_freq_rest(:, i), 'r');
    title(multiverse_labels{bb_pipeline_idx(i)});
end
ax = subplot(5, 3, 14);
ax.Visible = 'off';
lgd = legend(ax, [p1 p2], {'Within-hemisphere', 'Across-hemisphere'}, ...
    'Location', 'northoutside');
title(lgd, 'Connections');
sgtitle('abs(Coh)');


% abs(LagCoh)
abslagcoh_freq_rest = results_conn.abslagcoh_freq(bb_pipelines, :, :, :, :, :);
within_abslagcoh_freq_rest = reshape((abslagcoh_freq_rest(:, :, :, :, 1, 2) + ...
    abslagcoh_freq_rest(:, :, :, :, 3, 4)) / 2, n_bb_pipelines, [], n_freqs);
within_abslagcoh_freq_rest = squeeze(mean(within_abslagcoh_freq_rest(:, analyzed_sessions > 0, :), 2))';
across_abslagcoh_freq_rest = reshape((abslagcoh_freq_rest(:, :, :, :, 1, 3) + ...
    abslagcoh_freq_rest(:, :, :, :, 1, 4) + abslagcoh_freq_rest(:, :, :, :, 2, 3) + ...
    abslagcoh_freq_rest(:, :, :, :, 2, 4)) / 4, n_bb_pipelines, [], n_freqs);
across_abslagcoh_freq_rest = squeeze(mean(across_abslagcoh_freq_rest(:, analyzed_sessions > 0, :), 2))';

figure;
for i = 1:n_bb_pipelines
    subplot(5, 3, i); hold on;
    p1 = plot(freqs, within_abslagcoh_freq_rest(:, i), 'b');
    p2 = plot(freqs, across_abslagcoh_freq_rest(:, i), 'r');
    title(multiverse_labels{bb_pipeline_idx(i)});
end
ax = subplot(5, 3, 14);
ax.Visible = 'off';
lgd = legend(ax, [p1 p2], {'Within-hemisphere', 'Across-hemisphere'}, ...
    'Location', 'northoutside');
title(lgd, 'Connections');
sgtitle('abs(LagCoh)');

%% Draft figure for the manuscript
figure;
% pipelines_to_plot = 1:n_bb_pipelines;
pipelines_to_plot = [1 2 3];
colors = {'b', 'b', 'b', 'c', 'c', 'c', 'r', 'r', 'r', 'm', 'm', 'm'};
styles = {'-', '--', ':', '-', '--', ':', '-', '--', ':', '-', '--', ':'};
subplot(2, 3, 1); hold on;
for i = pipelines_to_plot    
    plot(freqs, within_absicoh_freq_rest(:, i), 'Color', colors{i}, 'LineStyle', styles{i});
end
xlim([3 40]); xlabel('Frequency (Hz)'); ylabel('abs(ImCoh)');
subplot(2, 3, 2); hold on;
for i = pipelines_to_plot  
    plot(freqs, within_abslagcoh_freq_rest(:, i), 'Color', colors{i}, 'LineStyle', styles{i});
end
xlim([3 40]); xlabel('Frequency (Hz)'); ylabel('abs(LagCoh)');
subplot(2, 3, 3); hold on;
for i = pipelines_to_plot  
    plot(freqs, within_coh_freq_rest(:, i), 'Color', colors{i}, 'LineStyle', styles{i});
end
xlim([3 40]); xlabel('Frequency (Hz)'); ylabel('abs(Coh)');
subplot(2, 3, 4); hold on;
for i = pipelines_to_plot   
    plot(freqs, across_absicoh_freq_rest(:, i), 'Color', colors{i}, 'LineStyle', styles{i});
end
xlim([3 40]); xlabel('Frequency (Hz)'); ylabel('abs(ImCoh)');
subplot(2, 3, 5); hold on;
for i = pipelines_to_plot   
    plot(freqs, across_abslagcoh_freq_rest(:, i), 'Color', colors{i}, 'LineStyle', styles{i});
end
xlim([3 40]); xlabel('Frequency (Hz)'); ylabel('abs(LagCoh)');
subplot(2, 3, 6); hold on;
for i = pipelines_to_plot    
    plot(freqs, across_coh_freq_rest(:, i), 'Color', colors{i}, 'LineStyle', styles{i});
end
xlim([3 40]); xlabel('Frequency (Hz)'); ylabel('abs(Coh)');

%% Export to R in the long format
pipeline_idx = repmat(bb_pipeline_idx, 1, n_freqs)';
conn_data = [
    repmat(freqs, n_bb_pipelines, 1), within_absicoh_freq_rest(:), ones(n_bb_pipelines * n_freqs, 1), ones(n_bb_pipelines * n_freqs, 1), pipeline_idx(:), multiverse_data(pipeline_idx(:), :); ...
    repmat(freqs, n_bb_pipelines, 1), within_abslagcoh_freq_rest(:), 2 * ones(n_bb_pipelines * n_freqs, 1), ones(n_bb_pipelines * n_freqs, 1), pipeline_idx(:), multiverse_data(pipeline_idx(:), :); ...
    repmat(freqs, n_bb_pipelines, 1), within_coh_freq_rest(:), 3 * ones(n_bb_pipelines * n_freqs, 1), ones(n_bb_pipelines * n_freqs, 1), pipeline_idx(:), multiverse_data(pipeline_idx(:), :); ...
    repmat(freqs, n_bb_pipelines, 1), across_absicoh_freq_rest(:), ones(n_bb_pipelines * n_freqs, 1), 2 * ones(n_bb_pipelines * n_freqs, 1), pipeline_idx(:), multiverse_data(pipeline_idx(:), :); ...
    repmat(freqs, n_bb_pipelines, 1), across_abslagcoh_freq_rest(:), 2 * ones(n_bb_pipelines * n_freqs, 1), 2 * ones(n_bb_pipelines * n_freqs, 1), pipeline_idx(:), multiverse_data(pipeline_idx(:), :); ...
    repmat(freqs, n_bb_pipelines, 1), across_coh_freq_rest(:), 3 * ones(n_bb_pipelines * n_freqs, 1), 2 * ones(n_bb_pipelines * n_freqs, 1), pipeline_idx(:), multiverse_data(pipeline_idx(:), :) ...
];
conn_labels = {'Freqs', 'Conn', 'Measure', 'Type', 'Pipeline', 'Inverse', 'Band', 'ROI_Agg', 'Mask', 'ROI_Method'};
pipeline_labels = multiverse_labels(bb_pipeline_idx);

%% Multiverse - Clear to Save RAM
clear results_conn;

%% Export the results to R
save([savedata 'BCI_MI_sanity_checks.mat'], ...
    'demo_snr', 'example_snr', 'mu_power_diff', ...
    'conn_data', 'conn_labels', 'pipeline_labels');