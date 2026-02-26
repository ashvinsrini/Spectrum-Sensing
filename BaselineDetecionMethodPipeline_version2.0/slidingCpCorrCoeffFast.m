function rhoVec = slidingCpCorrCoeffFast(x, Fs, L, Nwin, hop)
% normalized magnitude correlation |sum x1*conj(x2)| / sqrt(E1 E2)
    x = x(:);
    N = numel(x);
    starts = 1:hop:(N - (Nwin + L) + 1);
    rhoVec = zeros(numel(starts),1);
    for k = 1:numel(starts)
        s = starts(k);
        x1 = x(s:s+Nwin-1);
        x2 = x(s+L:s+L+Nwin-1);
        num = sum(x1 .* conj(x2));
        den = sqrt(sum(abs(x1).^2) * sum(abs(x2).^2) + eps);
        rhoVec(k) = abs(num) / den;
    end
end