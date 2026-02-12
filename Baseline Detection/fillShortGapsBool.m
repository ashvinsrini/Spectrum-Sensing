function y = fillShortGapsBool(x, maxGap)
% Fill short 0-gaps between 1-runs (boolean morphology).
    y = x(:) > 0;
    if maxGap <= 0 || isempty(y), return; end
    n = numel(y); i = 1;
    while i <= n
        if y(i) ~= 0
            i = i + 1; continue;
        end
        j = i;
        while j <= n && y(j) == 0
            j = j + 1;
        end
        gapLen = j - i;
        left  = (i>1) && y(i-1);
        right = (j<=n) && y(j);
        if gapLen <= maxGap && left && right
            y(i:j-1) = true;
        end
        i = j;
    end
end
