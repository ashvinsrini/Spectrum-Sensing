function A = readMaskAny(p)
% Reads uint8 mask from .hdf or .png (or other imread-supported formats)
    [A,~] = imread(p);
    if ndims(A) == 3
        error('Mask appears RGB at %s. Use the raw label mask (not maskenhanced).', p);
    end
    A = uint8(A);
end