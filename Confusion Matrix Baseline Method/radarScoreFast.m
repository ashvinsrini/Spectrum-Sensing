
function scoreRadar = radarScoreFast(x, Fs, PRI_range_s, PW_range_s, FsEnv)
% envelope downsample then call your radarPulseTrainScore
    x = x(:);
    D  = max(1, round(Fs/FsEnv));
    FsD = Fs / D;
    e = abs(x).^2;
    eSm = movmean(e, D);
    eD  = eSm(1:D:end);
    xEnv = sqrt(eD);
    [scoreRadar, ~, ~, ~] = radarPulseTrainScore(xEnv, FsD, PRI_range_s, PW_range_s);
end
