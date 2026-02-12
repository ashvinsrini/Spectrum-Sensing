function xw = applyTaper(x, Nfade)
% Raised-cosine fade in/out to reduce splatter at segment boundaries.
    x = x(:);
    N = numel(x);
    Nfade = min(Nfade, floor(N/2));
    if Nfade <= 1
        xw = x; return;
    end
    w = ones(N,1,'like',x);
    n = (0:Nfade-1).';
    ramp = 0.5*(1 - cos(pi*(n/(Nfade-1))));     % 0 -> 1
    w(1:Nfade) = ramp;
    w(end-Nfade+1:end) = flipud(ramp);
    xw = x .* w;
end