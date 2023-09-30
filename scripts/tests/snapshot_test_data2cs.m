rng('default');

%% Simulate some alpha oscillations mixed with noise
pnts = 1000;
nep = 100;
nchan = 5;
fs = 250;
[b, a] = butter(2, 2 * [8 12] / fs);
nfft = pnts;
n_overlap = fs;
peakbin = 13:19; % around 10 Hz

data = zeros(nchan, pnts, nep);
data_nb = zeros(nchan, pnts, nep);
for j = 1:nep
    y = filtfilt(b, a, rand(pnts, 1));
    [spec_sig, f] = pwelch(y', nfft, n_overlap, nfft, fs);
    
    fr_band_bins = dsearchn(f, [8; 12]);
    p_sig = sum(spec_sig(fr_band_bins(1):fr_band_bins(2), :));
    amp_sig = sqrt(1 ./ p_sig);
    y = y .* amp_sig';
    
    tot = zeros(nchan, pnts);
    for i = 1:5
        tot(i, :) = circshift(y, 3*(i-1))';
    end

    data(:, :, j) = tot;
end

%% Ground Truth Values
[ci, ri] = meshgrid(1:5);
cohy_gt = exp(3i * (ci - ri) * 10 * 2 * pi / fs);

%% Broadband
cs_bb = data2cs_bb(data, fres, 1);
cohy_bb = cs2coh(cs_bb);

%% Narrowband
cs_nb = data2cs_nb(data, 1);
cohy_nb = cs2coh(cs_nb);

%% Compare the results
cohs = {
    cohy_gt, ...
    squeeze(mean(cohy_bb(peakbin, :, :), 1)), ...
    squeeze(cohy_nb)
};
labels = {'gt', 'bb', 'nb'};
figure;
for k = 1:3
    subplot(1, 3, k);
    
    coh_vals = cohs{k};
    for i = 1:nchan
        for j = 1:nchan
            if (i == j)
                continue
            end
            
            con_value = coh_vals(i, j);
            polarplot([0 con_value]); hold on;
            text(angle(con_value), abs(con_value), [num2str(i) num2str(j)], ...
                'HorizontalAlignment','center',...
                'VerticalAlignment','bottom');
        end
    end
    
    rlim([0 1]);
    title(labels{k});
end
pause(1);

reply = input('Do you want to update the snapshot? y/n: ', 's');
switch (reply)
    case 'y'
        save([mfilename('fullpath') '.mat'], 'data', 'fres', 'cs_bb', 'cs_nb');
        fprintf('Snapshot was updated.\n');
    case 'n'
        fprintf('Snapshot was not updated.\n');
    otherwise
        fprintf('Snapshot was not updated. Unknown input.\n');
end 