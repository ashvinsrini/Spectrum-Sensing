
function out = mergeCloseBands(bandsHz, gapHz)
% Merge adjacent bands separated by <= gapHz.
    if isempty(bandsHz), out = bandsHz; return; end
    b = sortrows(bandsHz,1);
    out = b(1,:);
    for k = 2:size(b,1)
        if b(k,1) - out(end,2) <= gapHz
            out(end,2) = max(out(end,2), b(k,2));
        else
            out = [out; b(k,:)]; %#ok<AGROW>
        end
    end
end
