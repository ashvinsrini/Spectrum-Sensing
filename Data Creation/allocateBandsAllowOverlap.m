function [f0Map, bandsHz, labelsTrue] = allocateBandsAllowOverlap(presentNames, BW, fMin, fMax)
% Allow OVERLAPPING center frequencies (no collision/guard checks).
% presentNames: string array like ["NR","LTE",...]
% BW: struct with fields BW.NR, BW.LTE, BW.WLAN, BW.RADAR (Hz)
% fMin,fMax: placement limits for CENTER +/- BW/2 (Hz)

K = numel(presentNames);

bandsHz    = zeros(K,2);
labelsTrue = strings(K,1);
f0Map      = struct();

for k = 1:K
    tech = char(presentNames(k));
    bw   = BW.(tech);

    % Choose f0 uniformly at random (independent across techs => overlaps possible)
    f0lo = fMin + bw/2;
    f0hi = fMax - bw/2;
    if f0hi <= f0lo
        error('Band too wide for the [fMin,fMax] range: %s', tech);
    end

    f0 = f0lo + (f0hi - f0lo)*rand;

    f0Map.(tech) = f0;
    bandsHz(k,:) = [f0 - bw/2, f0 + bw/2];
    labelsTrue(k)= string(tech);
end
end
