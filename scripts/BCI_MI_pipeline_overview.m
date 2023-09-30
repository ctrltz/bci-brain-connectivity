%%%
% Prepare panels for the pipeline overview figure
%
% Outputs:
% 1. fig2c-roi-weights.png
%%%

% Plot AVG-flip and SVD weights for one participant and one pipeline
i_subject = 9;
i_session = 1;
i_pipeline = 2; % weights were saved only for 3SVD

roi_weights_avgflip = zeros(n_voxels, 1);
roi_weights_svd = zeros(n_voxels, 1);
for i_roi = 1:numel(ROI_inds)
    roi_ind = ROI_inds(i_roi);
    
    voxels_roi = get_voxels_roi(sa, roi_ind, []);

    svd_weights = results_spec.w_rois{i_pipeline, i_subject, i_session, i_roi};

    % reverse the polarity for display purposes for easier comparison
    % between AVG-flip and SVD
    if signflip(voxels_roi)' * svd_weights(:, 1) < 0
        svd_weights(:, 1) = -svd_weights(:, 1);
    end

    roi_weights_avgflip(voxels_roi) = signflip(voxels_roi);
    roi_weights_svd(voxels_roi) = svd_weights(:, 1);
end

% Scale SVD weights to [-1, 1] range
roi_weights_svd = roi_weights_svd ./ max(abs(roi_weights_svd));

h = figure('Position', [10 10 600 500]);
% AVG-flip
ax = subplot(121);
allplots_cortex_subplots(sa, cort5K2full(roi_weights_avgflip, sa), [-1 1], ...
    cm17, 'CSP-L', 1, 'views', 5, 'ax', ax); 
text(0.5, 1.2, 'AVG-flip', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', 18, ...
    'FontName', 'Arial');
% CSP-masked ROI definitions
ax = subplot(122);
allplots_cortex_subplots(sa, cort5K2full(roi_weights_svd, sa), [-1 1], ...
    cm17, 'CSP-L', 1, 'views', 5, 'ax', ax); 
text(0.5, 1.2, '1SVD', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', 18, ...
    'FontName', 'Arial');
exportgraphics(h, [cfg.asset.base 'fig2c-roi-weights.png']);
if (run_cmd)
    close(h);
end

%% Connectivity Metrics
[X, Y] = meshgrid(-1:0.01:1, -1:0.01:1);
D = X .^ 2 + Y .^ 2;

C = D;
LC = abs(Y ./ sqrt(1 - X .^ 2));
IC = abs(Y);

C(D > 1.05) = NaN;
LC(D > 1.05) = NaN;
IC(D > 1.05) = NaN;

T = (0:1:360) / 180 * pi;
Xc = cos(T);
Yc = sin(T);

h = figure('Position', [10 10 1250 300]);
subplot(141); hold on;
surf(X, Y, IC, IC, 'EdgeColor', 'none', 'FaceColor', 'interp'); view(2);
plot3(Xc, Yc, ones(size(T)), 'k', 'LineWidth', 2);
xlabel('Re'); ylabel('Im'); clim([0 1]); axis square;
title('ImCoh', 'FontSize', 18, 'FontName', 'Roboto', 'FontWeight', 'normal');
subplot(142); hold on;
surf(X, Y, LC, LC, 'EdgeColor', 'none', 'FaceColor', 'interp'); view(2);
plot3(Xc, Yc, ones(size(T)), 'k', 'LineWidth', 2);
xlabel('Re'); ylabel('Im'); clim([0 1]); axis square;
title('LagCoh', 'FontSize', 18, 'FontName', 'Roboto', 'FontWeight', 'normal');
subplot(143); hold on;
surf(X, Y, C, C, 'EdgeColor', 'none', 'FaceColor', 'interp'); view(2);
plot3(Xc, Yc, ones(size(T)), 'k', 'LineWidth', 2);
xlabel('Re'); ylabel('Im'); clim([0 1]); axis square;
title('Coherence', 'FontSize', 18, 'FontName', 'Roboto', 'FontWeight', 'normal');
colormap(cm17a);
ax = subplot(144);
allplots_cortex_subplots(sa, cort5K2full(roi_weights_svd, sa), [-1 1], ...
    cm17, 'arb. units', 1, 'views', 9, 'ax', ax, ...
    'cb_fontsize', 16, 'cb_fontname', 'Roboto');
exportgraphics(h, [cfg.results.misc 'connectivity-metrics.png'], ...
    'Resolution', 500);
if (run_cmd)
    close(h);
end

%% Template for Within / Across Edges
h = figure('Position', [10 10 600 400]);
allplots_cortex_subplots(sa, cort5K2full(zeros(n_voxels, 1), sa), [-1 1], ...
    cm17, 'CSP-L', 1, 'views', 5);
exportgraphics(h, [cfg.results.misc 'cortex-dorsal-view.png'], ...
    'Resolution', 500);
if (run_cmd)
    close(h);
end