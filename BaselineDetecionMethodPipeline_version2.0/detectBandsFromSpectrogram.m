function bandsHz = detectBandsFromSpectrogram(Sdb, fHz)
% Robust stripe detection from time-averaged spectrum (returns [K x 2] bands in Hz).
    fHz = fHz(:);
    df = abs(fHz(2)-fHz(1));

    % average PSD across time
    P = mean(10.^(Sdb/10), 2);

    % smooth
    smoothBins = max(3, round(2e6/(df+eps)));
    Psm = movmean(P, smoothBins);

    % robust noise floor via low-percentile bins
    idxNoise = Psm <= prctile(Psm, 40);
    muN  = median(Psm(idxNoise));
    sigN = 1.4826 * mad(Psm(idxNoise), 1);   % robust scale

    % threshold
    kOcc = 8;   % stronger => fewer false stripes
    thr = muN + kOcc*sigN;
    occ = (Psm >= thr);

    % morphology
    minBandHz = 8e6;
    gapHz     = 3e6;
    minBins = max(1, round(minBandHz/df));
    gapBins = max(1, round(gapHz/df));

    occ = fillShortGapsBool(occ, gapBins);
    occ = removeShortRunsBool(occ, minBins);

    bandsHz = binsToBands(fHz, occ);

    % merge very close bands (prevents 1 stripe splitting into 2)
    mergeHz = 5e6;
    bandsHz = mergeCloseBands(bandsHz, mergeHz);
end



