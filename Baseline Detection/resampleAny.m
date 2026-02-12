
function y = resampleAny(x, FsIn, FsOut)
% Uses polyphase resample if available, else linear interpolation.
    x = x(:);
    if abs(FsIn - FsOut) < 1e-9
        y = x;
        return;
    end
    if exist('resample','file') == 2
        % rational approximation for FsOut/FsIn
        [p,q] = rat(FsOut/FsIn, 1e-12);
        y = resample(x, p, q);
    else
        % fallback: linear interp
        tIn  = (0:numel(x)-1).'/FsIn;
        tOut = (0:round(numel(x)*FsOut/FsIn)-1).'/FsOut;
        y = interp1(tIn, x, tOut, 'linear', 0);
    end
end
