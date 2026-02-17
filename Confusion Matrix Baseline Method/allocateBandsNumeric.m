function [f0Map, bandsTrueHz, labelsTrue] = allocateBandsNumeric(presentNames, BW, fMin, fMax, guardHz)
% Random non-overlapping center frequencies.
    used = zeros(0,2);
    f0Map = struct();
    bandsTrueHz = zeros(0,2);
    labelsTrue  = strings(0,1);

    order = presentNames(randperm(numel(presentNames)));

    for k = 1:numel(order)
        tech = order(k);
        bw = BW.(char(tech));
        tries = 0;

        while true
            tries = tries + 1;
            if tries > 500
                error('Band allocation failed: reduce BW/guard or widen [fMin,fMax].');
            end

            fc = (fMin + bw/2 + guardHz) + rand() * ((fMax - bw/2 - guardHz) - (fMin + bw/2 + guardHz));
            fL = fc - bw/2;
            fH = fc + bw/2;

            ok = true;
            for u = 1:size(used,1)
                if ~(fH < used(u,1)-guardHz || fL > used(u,2)+guardHz)
                    ok = false; break;
                end
            end

            if ok
                used(end+1,:) = [fL fH]; %#ok<AGROW>
                f0Map.(char(tech)) = fc;
                bandsTrueHz(end+1,:) = [fL fH]; %#ok<AGROW>
                labelsTrue(end+1,1)  = string(tech); %#ok<AGROW>
                break;
            end
        end
    end
end