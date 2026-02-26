function y = removeShortRunsBool(x, minLen)
    y = x(:) > 0;
    if minLen <= 1 || isempty(y), return; end
    n = numel(y); i = 1;
    while i <= n
        val = y(i);
        j = i;
        while j <= n && y(j) == val, j = j + 1; end
        runLen = j - i;
        if val == 1 && runLen < minLen
            y(i:j-1) = false;
        end
        i = j;
    end
end