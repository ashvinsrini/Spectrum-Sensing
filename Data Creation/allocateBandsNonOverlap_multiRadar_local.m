function [labelsTrue, bandsTrueHz, f0ListHz] = allocateBandsNonOverlap_multiRadar_local( ...
    presentOther, BW, fMin, fMax, guardHz, numRadarRegions)
% Allocate non-overlapping frequency bands for:
%   presentOther (subset of NR/LTE/WLAN) + RADAR repeated numRadarRegions times
% Ensures:
%   - RADAR regions do not overlap each other
%   - RADAR regions do not overlap any other tech
%   - all separations have at least guardHz margin (implemented as +/- guardHz/2)

labelsTrue = [presentOther(:); repmat("RADAR", numRadarRegions, 1)];
K = numel(labelsTrue);

bwList = zeros(K,1);
for k = 1:K
    bwList(k) = BW.(labelsTrue(k));
end

% Place wider bands first for higher success probability
effW = bwList + guardHz;
[~,ord] = sort(effW, 'descend');

maxRestarts = 200;
maxTriesPerBand = 2000;

for restart = 1:maxRestarts
    f0Tmp    = zeros(K,1);
    bandTmp  = zeros(K,2);
    occExp   = zeros(0,2);   % occupied expanded bands: [fL-guard/2, fH+guard/2]
    ok = true;

    for ii = 1:K
        k = ord(ii);
        bw = bwList(k);

        f0min = fMin + bw/2;
        f0max = fMax - bw/2;
        if f0max <= f0min
            ok = false; break;
        end

        placed = false;
        for tr = 1:maxTriesPerBand
            f0 = f0min + (f0max - f0min) * rand();
            band = [f0 - bw/2, f0 + bw/2];
            expBand = [band(1) - guardHz/2, band(2) + guardHz/2];

            if isempty(occExp)
                placed = true;
            else
                overlap = any(~(expBand(2) < occExp(:,1) | expBand(1) > occExp(:,2)));
                placed = ~overlap;
            end

            if placed
                f0Tmp(k)   = f0;
                bandTmp(k,:) = band;
                occExp = [occExp; expBand]; 
                break;
            end
        end

        if ~placed
            ok = false; break;
        end
    end

    if ok
        f0ListHz   = f0Tmp;
        bandsTrueHz = bandTmp;
        return;
    end
end

error('allocateBandsNonOverlap_multiRadar_local: failed to place non-overlapping bands. Try reducing numRadarRegions or guardHz.');
end
