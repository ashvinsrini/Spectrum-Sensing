function xw = widenRadarPulses_local(x, Fs, widen_us)
%WIDENRADARPULSES_LOCAL  Make radar pulses thicker in time by dilating pulse support.
%   x        : complex radar template (Nx1)
%   Fs       : sample rate (Hz)
%   widen_us : dilation half-width in microseconds (extends each pulse by +/- widen_us)

    x  = x(:);
    if widen_us <= 0
        xw = x;
        return;
    end

    a = abs(x);
    if ~any(a)
        xw = x;
        return;
    end

    % Robust pulse threshold: pick "high-energy" samples as pulse core
    thr = max(prctile(a, 99), 5*median(a + eps));
    core = (a > thr);

    % Fallback if core is too sparse
    if nnz(core) < 200
        thr  = prctile(a, 99.5);
        core = (a > thr);
    end

    if nnz(core) < 50
        % If still too few, give up rather than corrupting the template
        xw = x;
        return;
    end

    L = max(1, round(widen_us*1e-6*Fs));   % samples to extend on each side

    % Fast 1-D dilation using moving max (O(N))
    wide = movmax(uint8(core), [L L]) > 0;

    % Newly added samples
    add = wide & ~core;

    xw = x;

    % Fill added samples by recycling existing pulse samples (keeps band-consistent content)
    idxCore = find(core);
    idxAdd  = find(add);
    if ~isempty(idxAdd)
        pick = idxCore(randi(numel(idxCore), numel(idxAdd), 1));
        xw(idxAdd) = x(idxCore(pick));
    end
end
