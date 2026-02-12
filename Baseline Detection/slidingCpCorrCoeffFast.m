%% ==========================================================
% fast sliding-window normalized CP-correlation coefficient
%   rho[i] = |sum x[n] x*[n+L]|^2 / (sum|x[n]|^2 * sum|x[n+L]|^2)
%  Computed efficiently with cumulative sums (O(N))
%% ==========================================================
function [rho, tFeat] = slidingCpCorrCoeffFast(x, Fs, L, Nwin, hop)
    x = x(:);
    N = numel(x);
    if Nwin < 8, Nwin = 8; end
    if hop  < 1, hop  = 1; end
    if L <= 0 || (L + Nwin + 2) > N
        rho = zeros(0,1); tFeat = zeros(0,1); return;
    end

    % aligned sequences for lag L
    x1 = x(1:N-L);
    x2 = x(1+L:N);

    r  = x1 .* conj(x2);               % elementwise correlation sequence
    p1 = abs(x1).^2;
    p2 = abs(x2).^2;

    % sliding sums via cumulative sums
    cs_r  = cumsum([0; r]);
    cs_p1 = cumsum([0; p1]);
    cs_p2 = cumsum([0; p2]);

    % window starts i=1..M
    M = numel(r) - Nwin + 1;
    if M < 1
        rho = zeros(0,1); tFeat = zeros(0,1); return;
    end

    sum_r  = cs_r(1+Nwin:end)  - cs_r(1:end-Nwin);
    sum_p1 = cs_p1(1+Nwin:end) - cs_p1(1:end-Nwin);
    sum_p2 = cs_p2(1+Nwin:end) - cs_p2(1:end-Nwin);

    idx = 1:hop:M;                      % hop-sampled windows

    num = abs(sum_r(idx)).^2;
    den = (sum_p1(idx) .* sum_p2(idx)) + eps;

    rho = num ./ den;                   % normalized coefficient (0..~1)
    tFeat = ((idx - 1) + Nwin/2).' / Fs; % seconds (center of window)
end