function bandsHz = binsToBands(fHz, occ)
    occ = occ(:) > 0;
    fHz = fHz(:);
    n = numel(occ);
    bandsHz = zeros(0,2);
    i = 1;
    while i <= n
        if ~occ(i), i = i + 1; continue; end
        j = i;
        while j <= n && occ(j), j = j + 1; end
        bandsHz(end+1,:) = [fHz(i) fHz(j-1)]; %#ok<AGROW>
        i = j;
    end
end