function x = simpleLTElikeOFDM(Fs)
% Fallback LTE-like OFDM waveform (roughly 20 MHz, 15 kHz SCS)
% Not standards-accurate LTE, but OFDM-like with LTE-ish parameters.
    Nfft = 2048;
    Nsc  = 1200;                 % 100 PRB * 12
    Nsym = 14*10;                % ~10 subframes worth of OFDM symbols (arbitrary)
    Ncp  = 144;                  % LTE normal CP (typical first symbol differs; simplified)

    M = 16;
    data = qammod(randi([0 M-1], Nsc*Nsym, 1), M, 'UnitAveragePower', true);
    grid = reshape(data, [Nsc Nsym]);

    % Map to center of FFT bins
    X = zeros(Nfft, Nsym);
    k0 = floor(Nfft/2) - floor(Nsc/2) + 1;
    X(k0:k0+Nsc-1, :) = grid;

    % IFFT per symbol + CP
    xSym = ifft(ifftshift(X,1), Nfft, 1);
    xCP  = [xSym(end-Ncp+1:end,:); xSym];
    x    = xCP(:);

    % Resample to requested Fs (if needed)
    FsNom = 30.72e6;
    if abs(Fs - FsNom) > 1
        x = resampleAny(x, FsNom, Fs);
    end
end