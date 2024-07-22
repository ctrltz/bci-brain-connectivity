%%%
% Estimate SNR for all ROIs and pipelines in the multiverse analysis
%
% Outputs:
%
% 1. BCI_MI_SNR_rois_multiverse.mat - results_snr structure with the
% following fields:
%
%    methods - method for estimating SNR
%    fit_range - frequency range used to fit FOOOF
%    mu_band - frequency band for estimating SNR
%    snr_rois_mu [pipelines x subjects x sessions x ROIs] - estimated SNR
%    r2_rois_fooof [pipelines x subjects x sessions x ROIs] - r2 of FOOOF
%    snr_avg [pipelines x subjects x sessions] - SNR averaged over ROIs
%%%

%% Calculate the SNR (ROIs)
snr_rois_mu = NaN(n_multiverse, n_subjects, n_sessions, n_rois);
r2_rois_fooof = NaN(n_multiverse, n_subjects, n_sessions, n_rois);
fooof_fits = cell(n_multiverse, n_subjects, n_sessions);
snr_rois_orig = cell(n_multiverse, n_subjects, n_sessions);
parfor (subject = 1:n_subjects, num_workers)
% for subject = 1:3
    for session = 1:n_sessions
%     for session = 1
        if results_spec.analyzed_sessions(subject, session)
            fprintf('%d %d ', subject, session);
            for i_pipeline = 1:n_multiverse
                psd_roi = results_spec.spec_rois{i_pipeline, subject, session};
                n_chans = size(psd_roi, 1);
                snr_roi_mu = zeros(n_chans, 1);
                r2_fooof = zeros(n_chans, 1);
                fooof_roi = cell(n_chans, 1);

                for ch = 1:n_chans
                    % Run FOOOF
                    fooof_results = fooof(results_spec.freqs, squeeze(psd_roi(ch, :)), ...
                        fit_range, fooof_settings, return_model);
                    r2_fooof(ch) = fooof_results.r_squared;
                    fooof_roi{ch} = fooof_results;

                    snr_roi_mu(ch, :) = calculate_SNR_fooof(fooof_results, methods, mu_band);
                end

                snr_rois_mu(i_pipeline, subject, session, :) = mean(reshape(snr_roi_mu, n_rois, []), 2);
                r2_rois_fooof(i_pipeline, subject, session, :) = mean(reshape(r2_fooof, n_rois, []), 2);
                fooof_fits{i_pipeline, subject, session} = reshape(fooof_roi, n_rois, []);
                snr_rois_orig{i_pipeline, subject, session} = reshape(snr_roi_mu, n_rois, []);

                fprintf('.');
            end
            fprintf('\n');
        end
    end
end
fprintf('Done\n');

%% Average SNR over all ROIs
snr_avg = squeeze(mean(snr_rois_mu, 4)); 

%% Save the results
results_snr = struct();
results_snr.methods = methods;
results_snr.fit_range = fit_range;
results_snr.mu_band = mu_band;
results_snr.snr_rois_mu = snr_rois_mu;
results_snr.r2_rois_fooof = r2_rois_fooof;
results_snr.snr_avg = snr_avg;

save([savedata 'BCI_MI_SNR_rois_multiverse.mat'], 'results_snr', '-v7.3');
