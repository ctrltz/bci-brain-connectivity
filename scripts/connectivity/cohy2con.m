function [conn] = cohy2con(cohy, measure, use_abs)
% COHY2CON Calculate various coherency-based phase synchronization measures
%
% Parameters:
%   cohy - coherency matrix
%   measure - which measure to calculate
%       'icoh' - imaginary part of coherency
%       'coh' - coherence / absolute part of coherency 
%       'lagcoh' - lagged coherence
%       'angle' - phase lag (rectified to [-pi/2; pi/2])
%       'angle2' - phase lag (in [-pi; pi])
%   use_abs - if 1, return the absolute value of the result
%
% Results:
%   conn - phase synchronization estimates, same shape as cohy

    switch (measure)
        case 'icoh'
            conn = imag(cohy);
        case 'coh'
            conn = abs(cohy);
        case 'lagcoh'
            conn = imag(cohy) ./ sqrt(1 - real(cohy) .^ 2);
            % Fix NaN values due to zero division
            conn(isnan(conn) & ~isnan(cohy)) = 0;
        case 'angle'  % [-pi/2; pi/2]
            conn = atan(imag(cohy) ./ real(cohy));
        case 'angle2' % [-pi; pi]
            conn = angle(cohy);
        otherwise
            error('Bad measure in cohy2con');
    end

    if (use_abs)
        conn = abs(conn);
    end
end