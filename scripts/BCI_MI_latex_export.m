%%% 
% Export settings and CSP results for including them in the paper in TeX
%
% Outputs (stored in results/tex):
% 1. settings.tex - analysis parameters from BCI_MI_config.m
% 2. csp.tex - number of task-relevant sources for each ROI
%%%

total_sessions = sum(cellfun(@(x) ~isempty(x), preproc_info), [1 2]);
excluded_sessions = total_sessions - sum(good_sessions, [1 2]);

settings_to_export = {
    "% Dataset", [], []; ...
    "numSubjects", n_subjects, '%d'; ...
    "numFollowUpSessionsMin", n_sessions - 5, '%d'; ...
    "numFollowUpSessionsMax", n_sessions - 1, '%d'; ...

    "% Preprocessing", [], []; ...
    "numSessionsExcluded", excluded_sessions, "%d"; ...
    "downsampleFreq", downsample, '%d'; ...
    "mirrorSeconds", n_mirror_seconds, '%d'; ...

    "% Frequency band of interest", [], []; ... 
    "muLow", mu_band(1), '%d'; ...
    "muHigh", mu_band(2), '%d'; ...

    "% FOOOF", [], []; ...
    "fooofFitRangeLow", fit_range(1), '%d'; ...
    "fooofFitRangeHigh", fit_range(2), '%d'; ...
    "fooofPeakWidthMin", fooof_settings.peak_width_limits(1), '%d'; ...
    "fooofPeakWidthMax", fooof_settings.peak_width_limits(2), '%d'; ...
    "fooofNumPeaksMax", fooof_settings.max_n_peaks, '%d'; ...

    "% FFT", [], []; ...
    "fftFreqRes", srate / nfft, "%.2f"; ...
    "fftSegLength", nfft / srate, "%.1f"; ...

    "% Multiverse", [], []; ...
    "numBands", n_bands, '%d'; ...
    "numInverse", n_inverse, '%d'; ...
    "numROIDefs", n_roi_defs, '%d'; ...
    "numROIMethods", n_roi_methods, '%d'; ...
    "numPipelines", n_multiverse, '%d'; ...
    "numBBPipelines", n_multiverse / n_bands, '%d'; ...
    
    "% Inverse", [], []; ...
    "numVoxels", n_voxels, '%d'; ...
    "inverseRegFactor", lambda, '%.2f'; ...

    "% CSP", [], []; ...
    "cspSourceThreshold", csp_source_threshold, "%.1f"
};
n_settings = size(settings_to_export, 1);

fileID = fopen([cfg.results.tex 'settings.tex'], 'w');
for i = 1:n_settings
    name = settings_to_export{i, 1};
    value = settings_to_export{i, 2};
    format = settings_to_export{i, 3};

    if (isempty(value))
        fprintf(fileID, '\n%s\n', name);
    else
        fprintf(fileID, '%s\n', export_to_latex(name, value, format));
    end
end
fclose(fileID);

%% Export results to include them in LaTeX
% Report number of selected voxels per ROI sorted in descending order
% Highlight ROIs that were analyzed in bold
cspStr = "";
[~, sort_idx] = sort(num_voxels_roi, 'descend');
for i = 1:n_rois_HO
    if num_voxels_roi(sort_idx(i)) == 0
        break
    end

    if ismember(sort_idx(i), ROI_inds)
        newStr = sprintf('\\textbf{%d} & \\textbf{%s} \\\\ \n', ...
            num_voxels_roi(sort_idx(i)), sa.HO_labels{sort_idx(i)});
    else
        newStr = sprintf('%d & %s \\\\ \n', ...
            num_voxels_roi(sort_idx(i)), sa.HO_labels{sort_idx(i)});
    end
    cspStr = cspStr + newStr;
end

fileID = fopen([cfg.results.tex 'csp.tex'], 'w');
fprintf(fileID, '%s\n', export_to_latex('cspSourcesPerROI', cspStr, []));
fclose(fileID);
