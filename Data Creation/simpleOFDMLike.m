function x = simpleOFDMLike(Fs, Nfft, Nsc, Nsym)
% basic OFDM with CP (not standard LTE mapping; good enough baseline)
    M = 16;
    grid = qammod(randi([0 M-1], Nsc*Nsym, 1), M, 'UnitAveragePower', true);
    grid = reshape(grid, [Nsc Nsym]);

    % map to centered bins
    X = zeros(Nfft, Nsym);
    k0 = floor(Nfft/2) - floor(Nsc/2) + 1;
    X(k0:k0+Nsc-1,:) = grid;

    % CP lengths (rough LTE-like): first symbol longer, others shorter
    cp1 = 160;  cpr = 144;

    x = [];
    for s = 1:Nsym
        xt = ifft(ifftshift(X(:,s)), Nfft);
        cp = (s==1)*cp1 + (s>1)*cpr;
        x = [x; xt(end-cp+1:end); xt]; 
    end
end