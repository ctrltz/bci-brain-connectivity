%%% 
% Prepare the CSP figure and plot anatomical and task-based ROI
% definitions
%
% Outputs:
% 1. colors.mat - cm17 colormap and colors of ROI to use in R
% 2. fig3-group-csp-source-space-mask.png - CSP figure
% 3. fig2b-roi-definition.png - panel with anatomical and task-based ROIs
%%%

%% Prepare data for the figure
csp_rec = csp_a_group(sa.myinds, :)' * A_eloreta_normal;
lim = max(abs(csp_rec), [], 2);
prcthresh = prctile(abs(csp_rec), csp_source_threshold, 2);
mask = abs(csp_rec) > prcthresh;

assert(isequal(mask(1, :) | mask(end, :), ...
    csp_mask(thresholds == csp_source_threshold, :)), ...
    "Calculations of voxel mask do not match between CSP fit and plot");

% Scale CSP patterns and source reconstruction to [-1, 1] range
csp_a_disp = csp_a_group ./ max(abs(csp_a_group), [], 1);
csp_rec_disp = csp_rec ./ max(abs(csp_rec), [], 2);
mask_disp = mask(1, :) | mask(end, :);

% Make topomaps positive at C3/C4 for display
ind_C3 = find(ismember({all_chanlocs(:).labels}, {'C3'}));
if (csp_a_disp(ind_C3, 1) < 0)
    csp_a_disp(:, 1) = -csp_a_disp(:, 1);
end

ind_C4 = find(ismember({all_chanlocs(:).labels}, {'C4'}));
if (csp_a_disp(ind_C4, end) < 0)
    csp_a_disp(:, end) = -csp_a_disp(:, end);
end

% Calculate spectra limits in the displayed frequency range
freq_range = [2 35];
spec_csp_disp = spec_csp_avg(:, [1 end], ...
    freqs >= freq_range(1) & freqs <= freq_range(2));
hi = 10 * log10(max(spec_csp_disp, [], 'all'));
lo = 10 * log10(min(spec_csp_disp, [], 'all'));

%% Prepare the figure
myblue = cm17(1, :);
myred = cm17(end, :);
myalpha = 0.5;
h = figure('Position', [10, 10, 1030, 650]);
set(h, 'DefaultAxesFontName', 'Arial');
% CSP patterns
subplot(3, 4, 1);
topoplot(csp_a_disp(:, 1), all_chanlocs);
text(-0.1, 1, 'A', 'Units', 'normalized', ...
    'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
t = text(-0.25, 0.5, 'Right Hand', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'FontName', 'Arial');
t.Rotation = 90;
subplot(3, 4, 5);
topoplot(csp_a_disp(:, end), all_chanlocs);
t = text(-0.25, 0.5, 'Left Hand', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'FontName', 'Arial');
t.Rotation = 90;
% CSP spectra
subplot(3, 4, 2); hold on;
set(gca, 'FontSize', 10);
set(gca, 'layer', 'top');
set(gca, 'LineWidth', 1);
fill([mu_band(1) mu_band(1) mu_band(2) mu_band(2)], ...
    [lo hi hi lo], [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
p1 = plot(freqs, 10*log10(squeeze(spec_csp_avg(1, 1, :))), 'Color', myblue, 'LineWidth', 1);
p2 = plot(freqs, 10*log10(squeeze(spec_csp_avg(2, 1, :))), 'Color', myred, 'LineWidth', 1);
xlim(freq_range); ylim([lo hi]);
xticks([freq_range(1) 10 20 30 freq_range(2)]);
xlabel('Frequency (Hz)'); ylabel('10*log10(PSD)');
text(-0.3, 1, 'B', 'Units', 'normalized', ...
    'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
subplot(3, 4, 6); hold on;
set(gca, 'FontSize', 10);
set(gca, 'layer', 'top');
set(gca, 'LineWidth', 1);
fill([mu_band(1) mu_band(1) mu_band(2) mu_band(2)], ...
    [lo hi hi lo], [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
p1 = plot(freqs, 10*log10(squeeze(spec_csp_avg(1, end, :))), 'Color', myblue, 'LineWidth', 1);
p2 = plot(freqs, 10*log10(squeeze(spec_csp_avg(2, end, :))), 'Color', myred, 'LineWidth', 1);
xlim(freq_range); ylim([lo hi]);
xticks([freq_range(1) 10 20 30 freq_range(2)]);
xlabel('Frequency (Hz)'); ylabel('10*log10(PSD)');
ax = subplot(3, 4, 10);
ax.Visible = 'off';
ax.Position(2) = ax.Position(2) + ax.Position(4);
ax.Position(4) = 0;
lgd = legend(ax, [p1 p2], {'Right-hand', 'Left-hand'}, ...
    'Location', 'southoutside', 'FontSize', 9);
title(lgd, 'Imaginary Movement');
% Source reconstructed CSP patterns
ax = subplot(3, 4, 3);
allplots_cortex_subplots(sa, cort5K2full(abs(csp_rec_disp(1, :)), sa), [0 1], ...
    cm17a, 'CSP-R', 1, 'views', 5, 'ax', ax);
text(-0.3, 1, 'C', 'Units', 'normalized', ...
    'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
ax = subplot(3, 4, 7);
allplots_cortex_subplots(sa, cort5K2full(abs(csp_rec_disp(end, :)), sa), [0 1], ...
    cm17a, 'CSP-L', 1, 'views', 5, 'ax', ax);
% Thresholded source reconstructed CSP patterns
ax = subplot(3, 4, 4);
allplots_cortex_subplots(sa, cort5K2full(mask(1, :), sa), [0 1], ...
    cm17a, 'CSP-R', 1, 'views', 5, 'ax', ax);
text(-0.3, 1, 'D', 'Units', 'normalized', ...
    'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
ax = subplot(3, 4, 8);
allplots_cortex_subplots(sa, cort5K2full(mask(end, :), sa), [0 1], ...
    cm17a, 'CSP-L', 1, 'views', 5, 'ax', ax);
ax = subplot(3, 4, 11);
allplots_cortex_subplots(sa, cort5K2full(mask(end, :), sa), [-1 1], ...
    cm17, 'arb. units', 1, 'views', 9, 'ax', ax, 'cb_location', 'southoutside');
ax.Position(2) = ax.Position(2) + ax.Position(4);
ax.Position(4) = 0;
ax.Colorbar.Position = [0.565 0.275 0.11 0.02];
colormap(cm17); % re-paint CSP patterns
exportgraphics(h, [savepath_group_csp 'fig3-group-csp-source-space-mask.png'], ...
    'Resolution', 500);

%% Plots for the Figure 3
% Highlight voxels from the target ROIs with values from -1 to 1
roi_anatomical_mask = zeros(n_voxels, 1);
roi_csp_mask = zeros(n_voxels, 1);
values = [-0.5, -1, 0.5, 1];
pos = 1 + round((size(cm17, 1) - 1) * (values + 1) / 2);
for i_roi = 1:numel(ROI_inds)
    roi_ind = ROI_inds(i_roi);
    
    voxels_roi = get_voxels_roi(sa, roi_ind, []);
    mask_roi = get_voxels_roi(sa, roi_ind, csp_voxel_mask);

    roi_anatomical_mask(voxels_roi) = values(i_roi);
    roi_csp_mask(mask_roi) = values(i_roi);
end

h = figure('Position', [10 10 600 500]);
% Anatomical ROI definitions
ax = subplot(121);
allplots_cortex_subplots(sa, cort5K2full(roi_anatomical_mask, sa), [-1 1], ...
    cm17, 'CSP-L', 1, 'views', 5, 'ax', ax);
text(0.5, 1.2, 'Anatomical', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', 18, ...
    'FontName', 'Arial');
% CSP-masked ROI definitions
ax = subplot(122);
allplots_cortex_subplots(sa, cort5K2full(roi_csp_mask, sa), [-1 1], ...
    cm17, 'CSP-L', 1, 'views', 5, 'ax', ax);
text(0.5, 1.2, 'Task-based', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', 18, ...
    'FontName', 'Arial');
exportgraphics(h, [cfg.asset.base 'fig2b-roi-definition.png']);

%% Export the colors to use them in R
colors = struct();
colors.cm = cm17;
colors.pos = pos;

save([savedata 'colors.mat'], 'colors');