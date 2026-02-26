function [Sdb, f, t] = mySpectrogramDB(x, Fs, winLen, hop, nfft)
% Same as your earlier style, but returns:
%   Sdb: (Nf x Nt), f: (Nf x 1), t: (1 x Nt)
    x = x(:);
    w = hann(winLen,'periodic');

    N = numel(x);
    nFrames = floor((N - winLen)/hop) + 1;
    S = zeros(nfft, nFrames, 'single');
    t = zeros(1, nFrames);

    for m = 1:nFrames
        idx0 = (m-1)*hop + 1;
        seg  = x(idx0:idx0+winLen-1) .* w;
        X    = fftshift(fft(seg, nfft));
        S(:,m) = single(abs(X).^2);
        t(m) = (idx0 - 1 + winLen/2)/Fs;
    end

    f = (-nfft/2 : nfft/2-1).' * (Fs/nfft);
    Sdb = 10*log10(double(S) + 1e-12);
end