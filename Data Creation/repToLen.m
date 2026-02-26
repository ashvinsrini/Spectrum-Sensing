function y = repToLen(x, N)
    x = x(:);
    if numel(x) >= N
        y = x(1:N);
    else
        r = ceil(N/numel(x));
        y = repmat(x, r, 1);
        y = y(1:N);
    end
end