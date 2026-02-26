function x = radarPulseTrainBaseband(Fs, Ntotal)
% Simple LFM pulse train spanning the full 40 ms (baseband)
% Pulse width ~ [50,100] us, PRI ~ [1,2] ms, chirp BW ~ [5,20] MHz
    x = complex(zeros(Ntotal,1,'single'));
    T = Ntotal/Fs;

    PRI  = (1e-3 + rand()*(2e-3-1e-3));
    PW   = (50e-6 + rand()*(100e-6-50e-6));
    BWc  = (5e6 + rand()*(20e6-5e6));

    t0s = 0:PRI:(T-PRI);
    for k = 1:numel(t0s)
        t0 = t0s(k);
        n0 = round(t0*Fs)+1;
        L  = round(PW*Fs);
        if n0+L-1 > Ntotal, break; end

        tt = (0:L-1).'/Fs;
        mu = BWc/PW;   % chirp rate (Hz/s)
        % baseband LFM: exp(j*pi*mu*t^2)
        pulse = exp(1j*pi*mu*(tt.^2));
        x(n0:n0+L-1) = x(n0:n0+L-1) + single(pulse);
    end

    % mild taper to reduce splatter
    w = single(hann(min(Ntotal, 4096)));
    x(1:numel(w)) = x(1:numel(w)).*w;
end