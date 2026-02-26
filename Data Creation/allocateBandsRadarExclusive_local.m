function [labelsTrue, bandsTrueHz, f0ListHz] = allocateBandsRadarExclusive_local( ...
    presentOther, BW, fMin, fMax, guardHz, numRadarRegions, allowOverlapNonRadar)
% Places:
%  - numRadarRegions RADAR bands: non-overlapping with each other and exclusive vs all
%  - presentOther (NR/LTE/WLAN): if allowOverlapNonRadar=true, they may overlap each other
%    BUT they will never overlap any RADAR band.

edgeMarginHz = 12e6;  % helps avoid STFT weirdness at edges

% ---------- 1) Place RADAR bands non-overlapping ----------
labelsRadar = repmat("RADAR", numRadarRegions, 1);
Krad = numRadarRegions;

bwRadar = BW.RADAR;
occExp = zeros(0,2);    % expanded occupied bands: [fL-guard/2, fH+guard/2]

bandsRadar = zeros(Krad,2);
f0Radar    = zeros(Krad,1);

maxTries = 5000;
for k = 1:Krad
    placed = false;
    for tr = 1:maxTries
        f0min = (fMin + edgeMarginHz) + bwRadar/2;
        f0max = (fMax - edgeMarginHz) - bwRadar/2;
        f0 = f0min + (f0max - f0min)*rand();

        band = [f0 - bwRadar/2, f0 + bwRadar/2];
        expBand = [band(1) - guardHz/2, band(2) + guardHz/2];

        if isempty(occExp)
            ok = true;
        else
            ok = ~any(~(expBand(2) < occExp(:,1) | expBand(1) > occExp(:,2)));
        end

        if ok
            bandsRadar(k,:) = band;
            f0Radar(k)      = f0;
            occExp = [occExp; expBand]; 
            placed = true;
            break;
        end
    end
    if ~placed
        error('Failed to place %d non-overlapping RADAR bands. Reduce guardHz or numRadarRegions.', numRadarRegions);
    end
end

% ---------- 2) Place OTHER tech bands (overlap allowed among themselves) ----------
labelsOther = presentOther(:);
Koth = numel(labelsOther);

bandsOther = zeros(Koth,2);
f0Other    = zeros(Koth,1);

for k = 1:Koth
    tech = labelsOther(k);
    bw = BW.(tech);

    placed = false;
    for tr = 1:maxTries
        f0min = (fMin + edgeMarginHz) + bw/2;
        f0max = (fMax - edgeMarginHz) - bw/2;
        f0 = f0min + (f0max - f0min)*rand();

        band = [f0 - bw/2, f0 + bw/2];
        expBand = [band(1) - guardHz/2, band(2) + guardHz/2];

        % Must NOT overlap RADAR occupied regions
        okRadar = ~any(~(expBand(2) < occExp(:,1) | expBand(1) > occExp(:,2)));

        if ~okRadar
            continue;
        end

        if allowOverlapNonRadar
            % overlap among NR/LTE/WLAN is allowed -> accept
            placed = true;
        else
            % non-overlap among NR/LTE/WLAN also required:
            % add their exp bands into occExp temporarily to enforce
            placed = true;
            occExp = [occExp; expBand]; 
        end

        if placed
            bandsOther(k,:) = band;
            f0Other(k)      = f0;
            break;
        end
    end

    if ~placed
        error('Failed to place %s band without overlapping RADAR. Try smaller guardHz.', tech);
    end
end

% ---------- 3) Combine outputs ----------
labelsTrue  = [labelsOther; labelsRadar];
bandsTrueHz = [bandsOther; bandsRadar];
f0ListHz    = [f0Other; f0Radar];
end
