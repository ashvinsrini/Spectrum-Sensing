function [f0Map2,bands2,labels2] = enforceRadarNoOverlap_local( ...
    presentNames, BW, fMin, fMax, guardHz, f0Map, bandsTrueHz, labelsTrue)

f0Map2 = f0Map; bands2 = bandsTrueHz; labels2 = labelsTrue;

rid = find(labelsTrue=="RADAR",1);
if isempty(rid), return; end

rL = bandsTrueHz(rid,1); rH = bandsTrueHz(rid,2);

overlap = false;
for k = 1:size(bandsTrueHz,1)
    if k==rid, continue; end
    if bandsTrueHz(k,1) <= rH && bandsTrueHz(k,2) >= rL
        overlap = true;
        break;
    end
end

if overlap
    % Reallocate with non-overlap so radar is not hidden
    [f0Map2,bands2,labels2] = allocateBandsNumeric(presentNames, BW, fMin, fMax, guardHz);
end
end
