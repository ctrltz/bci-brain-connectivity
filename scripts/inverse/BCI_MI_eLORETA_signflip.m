%%% 
% Pre-compute eLORETA inverse operator and sign flips for averaging ROI
% time series
%
% Outputs:
% 1. precomputed/BCI_MI_sa_eLoreta.mat with the following variables:
%
%   clab - labels of all channels
%   sa - source analysis structure
%   L, L_normal - lead field (free and fixed orientation)
%   vc - voxel coordinates
%   signflip - vector of 1s and -1s reflecting whether time series of the
%       voxel should be flipped before averaging
%   all_chanlocs - chanlocs structure for all channels
%   n_sensors - number of channels
%   n_voxels - number of voxels
%   A_eloreta, A_eloreta_normal - inverse operator (free / fixed)
%   A_eloreta_normal_focal, A_eloreta_normal_smooth - inverse operator
%       (regularization parameter is different)
%%%

%% Working Directory
folder = 'precomputed/';

%% Dataset
dataset.bci_mi = 'BCI_MI';
dataset.lemon = 'LEMON';
dataset_to_use = dataset.bci_mi;
% dataset_to_use = dataset.lemon;

load cm17;

%% Prepare the leadfield matrix and inverse operator
if strcmp(dataset_to_use, dataset.bci_mi)
    load([folder 'BCI_MI_chanlocs60.mat']);
    all_chanlocs = chanlocs60;
    clab = {all_chanlocs(:).labels};
    % Fix case so that channel labels match the pre-computed ones
    clab([1, 2, 3, 10, 19, 28, 37, 46, 54, 59]) = ...
        {'Fp1', 'Fpz', 'Fp2', 'Fz', 'FCz', 'Cz', 'CPz', 'Pz', 'POz', 'Oz'};
    sa = prepare_sourceanalysis(clab, 'nyhead', struct('newlocs', 1));
elseif strcmp(dataset_to_use, dataset.lemon)
    load([folder 'LEMON_chanlocs61.mat']);
    all_chanlocs = chanlocs61;
    clab = {all_chanlocs(:).labels};
    sa = prepare_sourceanalysis(clab, 'nyhead', struct('newlocs', 1));
    sa.myinds = 1:numel(clab);
else
    error(['Unknown dataset: ' dataset_to_use]);
end

% common average reference transform
M = length(sa.clab_electrodes);
H = eye(M) - ones(M) / M;

% Lead field for 5K vertices restricted to the cortical surface with either
% free (L) or fixed along the normal to the cortex (L_normal) orientations
sa.voxels_5K_cort = sa.cortex5K.in_from_cortex75K(sa.cortex5K.in_cort);
L = sa.cortex75K.V_fem(:, sa.voxels_5K_cort, :);
L_normal = sa.cortex75K.V_fem_normal(:, sa.voxels_5K_cort);
[~, n_voxels, ~] = size(L);

% Apply common average reference to match preprocessing
L = reshape(H * reshape(L, M, []), M, n_voxels, 3);
L_normal = H * L_normal;

% Sizes
n_sensors = numel(all_chanlocs);
n_voxels = size(L, 2);

% Voxel coordinates
vc = sa.cortex75K.vc(sa.voxels_5K_cort, :);
vc_smooth = sa.cortex75K.vc_smooth(sa.voxels_5K_cort, :);

% eLoreta inverse operator
A_eloreta = mkfilt_eloreta2(L, lambda);
A_eloreta_normal = mkfilt_eloreta2(L_normal, lambda);
A_eloreta_normal_focal = mkfilt_eloreta2(L_normal, 0.001);
A_eloreta_normal_smooth = mkfilt_eloreta2(L_normal, 0.5);

% signflip mask per ROI
n_rois_HO = 96;  % 97 - 1 (subcortical)
signflip = zeros(numel(sa.voxels_5K_cort), 1);
for i_roi = 1:n_rois_HO
    voxels_roi = get_voxels_roi(sa, i_roi, []);
    
    normals_roi = sa.cortex75K.normals(sa.voxels_5K_cort(voxels_roi), :);

    % sign flip before average from brainstorm
    % Take the SVD to get the dominant orientation in this patch
    % v(:,1) is the dominant orientation
    [u, ~, ~] = svd(normals_roi, 0);
    % Get the flip mask for the data values
    flipmask = sign(u(:,1));
    % We want to flip the sign of the minimum number of time series, so
    % if there are mostly positive values: multiply the values by -FlipMask
    if (nnz(flipmask > 0) < nnz(flipmask < 0))
        flipmask = -flipmask;
    end
    
    signflip(voxels_roi) = flipmask;
end

% plot the signflip mask of the cortex surface
h = figure;
views = [5 9; 2 4];
allplots_cortex_subplots(sa, cort5K2full(signflip, sa), [-1 1], cm17, 'sign flip', 1, 'views', views);
exportgraphics(h, [cfg.results.misc 'signflip_cortex.png']);

save([folder dataset_to_use '_sa_eLoreta.mat'], ...
     'clab', 'sa', 'L', 'L_normal', 'vc', 'signflip', ...
     'all_chanlocs', 'n_sensors', 'n_voxels', ...
     'A_eloreta', 'A_eloreta_normal',  ...
     'A_eloreta_normal_focal', 'A_eloreta_normal_smooth');