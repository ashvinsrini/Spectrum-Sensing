
function [f0Map, bandsTrue] = allocateBands(presentNames, BW, fMin, fMax, guardHz)
% Random non-overlapping center frequencies for present technologies.
% Returns:
%   f0Map.(tech) = center freq
%   bandsTrue    = [K x 3] string table: [fL fH label]
    used = zeros(0,2);
    f0Map = struct();
    bandsTrue = strings(0,3);

    order = presentNames(randperm(numel(presentNames)));

    for k = 1:numel(order)
        tech = order(k);
        bw = BW.(tech);
        tries = 0;

        while true
            tries = tries + 1;
            if tries > 500
                error('Band allocation failed: increase guardHz or reduce BW.');
            end

            fc = (fMin + bw/2 + guardHz) + rand() * ((fMax - bw/2 - guardHz) - (fMin + bw/2 + guardHz));
            fL = fc - bw/2;
            fH = fc + bw/2;

            % check overlap with used intervals (with guard)
            ok = true;
            for u = 1:size(used,1)
                if ~(fH < used(u,1)-guardHz || fL > used(u,2)+guardHz)
                    ok = false; break;
                end
            end

            if ok
                used(end+1,:) = [fL fH]; 
                f0Map.(tech) = fc;

                bandsTrue(end+1,:) = [string(fL), string(fH), string(tech)]; 
                break;
            end
        end
    end

    % convert bandsTrue back to numeric bounds for plotting convenience later
    % keep it as strings but we parse at plot time
end