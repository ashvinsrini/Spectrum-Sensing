function [scoreRadar, priHat_s, duty, thr] = radarPulseTrainScore(x, Fs, PRI_range_s, PW_range_s)
% Detect pulsed radar by envelope periodicity in PRI_range_s.
% Returns:
%  scoreRadar in [0,1], estimated PRI (seconds), duty cycle, and threshold used.

    x = x(:);
    e = abs(x).^2;

    % Smooth envelope a bit (helps at low SNR)
    smooth_us = 10; % ~10 us
    Ns = max(1, round(smooth_us*1e-6*Fs));
    eS = movmean(e, Ns);

    % Robust threshold for "pulse present"
    mu  = median(eS);
    sig = 1.4826 * mad(eS, 1);     % robust scale
    kThr = 6;                      % tune 5..10
    thr  = mu + kThr*sig;

    b = (eS > thr);

    % Clean up mask: remove tiny runs, fill tiny gaps
    PWmin = PW_range_s(1); PWmax = PW_range_s(2);
    minRun = max(1, round(0.5*PWmin*Fs));   % keep >= half min pulse
    maxGap = max(1, round(0.2*PWmin*Fs));   % fill short gaps inside pulse

    b = fillShortGapsBool(b, maxGap);
    b = removeShortRunsBool(b, minRun);

    duty = mean(b);

    % If nothing detected, radar score = 0
    if nnz(b) < 10
        scoreRadar = 0; priHat_s = NaN; return;
    end

    % PRI search lags
    Lmin = max(1, round(PRI_range_s(1)*Fs));
    Lmax = min(numel(b)-2, round(PRI_range_s(2)*Fs));
    if Lmax <= Lmin
        scoreRadar = 0; priHat_s = NaN; return;
    end

    % Score(L) = overlap / (# active samples)
    denom = nnz(b);
    scores = zeros(Lmax-Lmin+1,1);
    Lgrid  = (Lmin:Lmax).';

    for ii = 1:numel(Lgrid)
        L = Lgrid(ii);
        scores(ii) = nnz( b(1:end-L) & b(1+L:end) ) / (denom + eps);
    end

    [scoreRadar, idx] = max(scores);
    priHat_s = Lgrid(idx)/Fs;
end
