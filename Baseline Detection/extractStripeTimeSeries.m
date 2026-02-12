%% ==========================================================
%   isolate stripe time series by mixing to baseband + lowpass
% ==========================================================
function [xBB_filt, fcHz] = extractStripeTimeSeries(x, Fs, fL, fH, firOrd)
    x = x(:);
    fcHz = 0.5*(fL + fH);
    bwHz = max(1, (fH - fL));    % avoid zero

    % Mix stripe to baseband
    n = (0:numel(x)-1).';
    xBB = x .* exp(-1i*2*pi*(fcHz/Fs)*n);

    % Lowpass filter to keep only this stripe (slight margin)
    cutoffHz = min(0.49*Fs, 0.6*(bwHz/2));               % 0.6*BW/2 keeps guard
    Wn = min(0.99, cutoffHz/(Fs/2));                     % normalized cutoff
    w  = hann(firOrd+1);
    h  = fir1(firOrd, Wn, 'low', w);

    % Filter (use filtfilt if available to avoid group delay)
    if exist('filtfilt','file') == 2
        xBB_filt = filtfilt(h, 1, xBB);
    else
        xBB_filt = filter(h, 1, xBB);
    end
end