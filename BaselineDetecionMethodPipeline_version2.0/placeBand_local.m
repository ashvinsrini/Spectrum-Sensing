function [fL,fH,f0] = placeBand_local(tech, BW, fMin, fMax, guardHz, forbiddenBands)
% Place one band for 'tech' uniformly, avoiding forbiddenBands (+guardHz).
% forbiddenBands: Kx2 [fL fH] (can be empty)

    BWk   = BW.(tech);
    f0min = fMin + BWk/2 + guardHz;
    f0max = fMax - BWk/2 - guardHz;

    if f0min >= f0max
        error('No room to place %s band. Reduce BW/guardHz or expand [fMin,fMax].', tech);
    end

    placed = false;
    for attempt = 1:800
        f0 = f0min + (f0max - f0min)*rand();
        fL = f0 - BWk/2;
        fH = f0 + BWk/2;

        if isempty(forbiddenBands)
            placed = true;
        else
            prev = forbiddenBands;
            placed = all( (fH < prev(:,1)-guardHz) | (fL > prev(:,2)+guardHz) );
        end

        if placed, break; end
    end

    if ~placed
        error('Could not place %s band without overlap. Reduce guardHz/BW or allowOverlapNonRadar.', tech);
    end
end
