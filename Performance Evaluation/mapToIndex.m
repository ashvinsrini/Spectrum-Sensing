function idx = mapToIndex(maskU8, labelVals)
% Converts mask values in labelVals -> indices 1..K
% Errors if any pixel has an unexpected value.
    K = numel(labelVals);
    idx = zeros(size(maskU8), 'uint16');

    for k = 1:K
        idx(maskU8 == labelVals(k)) = uint16(k);
    end

    if any(idx(:) == 0)
        bad = unique(maskU8(idx==0));
        error('Found unexpected label values in mask: %s', mat2str(double(bad(:)')));
    end
end