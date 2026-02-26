function bands2 = mergeCloseBands(bands, mergeHz)
    if isempty(bands), bands2 = bands; return; end
    bands = sortrows(bands,1);
    out = bands(1,:);
    for k = 2:size(bands,1)
        if bands(k,1) - out(end,2) <= mergeHz
            out(end,2) = max(out(end,2), bands(k,2));
        else
            out(end+1,:) = bands(k,:); 
        end
    end
    bands2 = out;
end