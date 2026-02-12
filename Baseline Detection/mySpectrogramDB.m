function [Sdb, f, t] = mySpectrogramDB(x, Fs, winLen, hop, nfft)
% STFT: returns Sdb as (freqBins x timeFrames)
    x = x(:);
    w = hannLocal(winLen);

    N = numel(x);
    nFrames = floor((N - winLen)/hop) + 1;
    if nFrames < 1
        error('Signal too short for spectrogram window.');
    end

    S = zeros(nfft, nFrames, 'single');
    t = zeros(1, nFrames);

    for m = 1:nFrames
        idx0 = (m-1)*hop + 1;
        seg  = x(idx0:idx0+winLen-1) .* w;
        X    = fftshift(fft(seg, nfft));
        S(:,m) = abs(X).^2;
        t(m) = (idx0 - 1 + winLen/2)/Fs;
    end

    f = (-nfft/2 : nfft/2-1) * (Fs/nfft);
    Sdb = 10*log10(S + 1e-12);
end

function w = hannLocal(N)
    n = (0:N-1).';
    w = 0.5 - 0.5*cos(2*pi*n/(N-1));
    w = single(w);
end