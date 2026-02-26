function xOut = forceRadarTwoStripes_local(xIn, Fs, nStripes)
% Keeps only nStripes high-energy pulse segments from an existing RADAR template.
% If fewer than nStripes segments exist, duplicates the strongest segment.

xIn = xIn(:);
N   = numel(xIn);

% Smooth envelope for robust pulse detection
env = abs(xIn);
Ls  = max(1, round(5e-6*Fs));   % 5 us smoothing
envS = movmean(env, Ls);

% Threshold: robust + also relative to peak
thr = max( median(envS) + 4*mad(envS,1), 0.25*max(envS) );

m = envS > thr;
d = diff([false; m; false]);
s = find(d==1);
e = find(d==-1)-1;

% Remove very short segments (too tiny to show as a stripe)
minLen = max(1, round(25e-6*Fs));  % 25 us
len = e - s + 1;
keep = len >= minLen;
s = s(keep); e = e(keep);

xOut = zeros(N,1,'like',xIn);

% If nothing detected, fall back to two fixed windows (still using xIn content)
if isempty(s)
    L  = max(1, round(120e-6*Fs));     % 120 us stripe
    gap = max(1, round(3e-3*Fs));      % 3 ms separation
    s1 = max(1, round(0.25*N));
    s2 = min(N-L+1, s1+gap);
    xOut(s1:s1+L-1) = xIn(s1:s1+L-1);
    xOut(s2:s2+L-1) = xIn(s2:s2+L-1);
    return;
end

% Segment energy
segE = zeros(numel(s),1);
for k = 1:numel(s)
    segE(k) = sum(abs(xIn(s(k):e(k))).^2);
end
[~,ord] = sort(segE,'descend');

% Pick segments with separation so we get distinct stripes
minGap = max(1, round(2e-3*Fs));  % 2 ms
sel = [];
for k = 1:numel(ord)
    kk = ord(k);
    if isempty(sel)
        sel(end+1) = kk; 
    else
        if all(abs(s(kk) - s(sel)) > minGap)
            sel(end+1) = kk; 
        end
    end
    if numel(sel) == nStripes
        break;
    end
end

% If still not enough, just take top ones
if numel(sel) < nStripes
    sel = ord(1:min(nStripes, numel(ord)));
end

% Copy selected segments
for k = 1:numel(sel)
    xOut(s(sel(k)):e(sel(k))) = xIn(s(sel(k)):e(sel(k)));
end

% If we only found 1 segment but need 2, duplicate it shifted in time
if numel(sel) < nStripes
    s1 = s(sel(1)); e1 = e(sel(1));
    L  = e1 - s1 + 1;
    gap = minGap;
    s2 = min(N-L+1, s1 + gap);
    xOut(s2:s2+L-1) = xIn(s1:e1);
end
end
