function y = repToLen(x, N)
% Repeat waveform x to length N and truncate.
    x = x(:);
    if numel(x) >= N
        y = x(1:N);
        return;
    end
    reps = ceil(N/numel(x));
    y = repmat(x, reps, 1);
    y = y(1:N);
end
