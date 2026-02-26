function maskNtNf = refineRadarMaskTF_local(maskNtNf, Sdb, trueOfFreq, labelVals)
% Robust RADAR mask refinement:
% - Process each RADAR vertical band separately
% - Detect pulse times from a robust time-trace (freq-percentile, not mean)
% - Adaptive threshold + FAILSAFE: keep top-K times if too few detected
% - Label full radar band at detected pulse times => horizontal stripes

noiseVal = uint8(labelVals(1));   % 0
radarVal = uint8(labelVals(5));   % 6

radarBinsAll = find(trueOfFreq == uint8(4));
if isempty(radarBinsAll), return; end

Nt = size(Sdb,2);

% ---- Split radar bins into contiguous segments (each segment = one RADAR vertical region) ----
d = diff(radarBinsAll);
cut = find(d > 1);
segStarts = [1; cut+1];
segEnds   = [cut; numel(radarBinsAll)];

% ---- Parameters (tune lightly if needed) ----
qFreq            = 90;    % percentile across freq inside radar band for time-trace (higher => more sensitive to pulses)
zThrInit         = 2.5;   % initial z-threshold on time-trace
minActiveFrac    = 0.01;  % at least 1% time bins should be active (else failsafe)
maxActiveFrac    = 0.40;  % avoid flooding
minActiveBinsAbs = 4;     % at least this many time bins active per RADAR band (failsafe)
topKFrac         = 0.03;  % failsafe: keep top ~3% time bins if detection too weak
timeDilateWin    = 5;     % thicken stripes in time
minRunLen        = 2;     % remove isolated 1-bin blips

for s = 1:numel(segStarts)
    seg = radarBinsAll(segStarts(s):segEnds(s));   % freq bins for this RADAR region

    % Reset entire RADAR region to Noise first
    maskNtNf(:, seg) = noiseVal;

    radarDb = double(Sdb(seg, :));                 % [nSeg x Nt]

    % ---- Robust time-trace: percentile over frequency bins (captures pulses better than mean) ----
    timeTrace = prctile(radarDb, qFreq, 1);        % [1 x Nt]

    % ---- Normalize: z-score using median + MAD (robust) ----
    mu  = median(timeTrace);
    sig = mad(timeTrace, 1) + 1e-9;
    z   = (timeTrace - mu) ./ sig;

    % ---- Adaptive threshold on z to keep active fraction in a sane range ----
    zThr = zThrInit;
    active = (z > zThr);
    for it = 1:10
        frac = mean(active);
        if frac < minActiveFrac
            zThr = zThr - 0.3;
        elseif frac > maxActiveFrac
            zThr = zThr + 0.3;
        else
            break;
        end
        active = (z > zThr);
    end

    % ---- FAILSAFE: if still too few active bins, keep the top-K strongest times ----
    nAct = nnz(active);
    minBins = max(minActiveBinsAbs, round(minActiveFrac * Nt));
    if nAct < minBins
        K = max(minBins, round(topKFrac * Nt));
        [~,ord] = sort(z, 'descend');
        active(:) = false;
        active(ord(1:min(K, Nt))) = true;
    end

    % ---- Thicken + clean (horizontal stripes) ----
    active = movmax(active, timeDilateWin) > 0;
    active = keepLongRuns1D_local(active, minRunLen);

    % Apply to mask (maskNtNf is Nt x Nf: time = rows)
    maskNtNf(active(:), seg) = radarVal;
end
end

function x = keepLongRuns1D_local(x, minLen)
% Keep only runs of true with length >= minLen (1D logical)
x = logical(x(:));
d = diff([false; x; false]);
s = find(d==1);
e = find(d==-1)-1;

x2 = false(size(x));
for k = 1:numel(s)
    if (e(k)-s(k)+1) >= minLen
        x2(s(k):e(k)) = true;
    end
end
x = x2;
end
