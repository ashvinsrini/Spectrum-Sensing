function y = resampleAny(x, FsIn, FsOut)
    x = x(:);
    if abs(FsIn - FsOut) < 1e-6
        y = x;
        return;
    end
    [p,q] = rat(FsOut/FsIn, 1e-10);
    y = resample(x, p, q);
end
