function maskNtNf = refineRadarMaskTF_slanted_local(maskNtNf, Sdb, trueOfFreq, labelVals, opts)
% Robust slanted-ridge RADAR refinement:
% - Uses "excess energy" relative to non-radar baseline
% - Time-gates to avoid always-on vertical fills
% - Enforces per-row TOP-K to avoid whole-band labeling
% - Auto-relaxes threshold if it detects nothing

    if nargin < 5 || isempty(opts), opts = struct; end

    if ~isfield(opts,'kMadStart'),     opts.kMadStart = 0.4; end
    if ~isfield(opts,'kActive'),       opts.kActive = 0.1; end
    if ~isfield(opts,'maxFracPerRow'), opts.maxFracPerRow = 0.08; end
    if ~isfield(opts,'minArea'),       opts.minArea = 6; end
    if ~isfield(opts,'ridgeLen'),      opts.ridgeLen = 7; end
    if ~isfield(opts,'dilateRad'),     opts.dilateRad = 1; end
    if ~isfield(opts,'useBothDirs'),   opts.useBothDirs = true; end

    radarId  = uint8(4);
    noiseVal = uint8(labelVals(1));
    radarVal = uint8(labelVals(5));

    Nf = numel(trueOfFreq);

    % ---- Make spectrogram as Nt x Nf ----
    if size(Sdb,1) == Nf
        Stf = Sdb.';   % (Nf x Nt) -> (Nt x Nf)
    else
        Stf = Sdb;     % already (Nt x Nf)
    end
    [Nt, Nf_check] = size(Stf);
    if Nf_check ~= Nf
        error('Mismatch: trueOfFreq=%d bins but spectrogram has %d bins.', Nf, Nf_check);
    end

    idxRadarF = (trueOfFreq(:) == radarId);
    if ~any(idxRadarF)
        return;
    end

    idxNonRadarF = ~idxRadarF;
    Srad = Stf(:, idxRadarF);           % Nt x Nr
    Nr   = size(Srad,2);

    % ---- Baseline per time from NON-radar frequencies (robust) ----
    if any(idxNonRadarF)
        Snoise = Stf(:, idxNonRadarF);
        base_t = median(Snoise, 2);     % Nt x 1
    else
        % fallback if only radar exists
        base_t = median(Srad, 2);
    end

    % "Excess energy" (dB) inside radar band relative to baseline
    E = bsxfun(@minus, Srad, base_t);   % Nt x Nr

    % small smoothing along frequency to reduce salt-pepper (optional but helps)
    E = movmedian(E, 3, 2);

    % ---- Time gating based on max excess per time ----
    maxE = max(E, [], 2);                       % Nt x 1
    medM = median(maxE);
    madM = mad(maxE, 1);                        % scalar (median abs deviation)
    sigM = 1.4826 * (madM + eps);
    thrActive = medM + opts.kActive * sigM;
    activeT = (maxE > thrActive);

    % If gate kills everything, keep top 20% most energetic times
    if ~any(activeT)
        p = prctile(maxE, 80);
        activeT = (maxE >= p);
    end

    % ---- Build BW with adaptive thresholds (auto-relax if empty) ----
    kList = [opts.kMadStart, opts.kMadStart-0.8, opts.kMadStart-1.6, 1.2];
    kList = kList(kList > 0.8);

    BWbest = false(Nt, Nr);
    bestScore = -inf;

    % Per-row top-K cap
    Kcap = max(1, ceil(opts.maxFracPerRow * Nr));
    Kcap = min(Kcap, 20);  % safety cap

    for kk = 1:numel(kList)
        kMad = kList(kk);

        % per-time robust threshold inside radar band using excess energy
        medE = median(E, 2);
        madE = mad(E, 1, 2);
        sigE = 1.4826 * (madE + eps);
        thrE = medE + kMad * sigE;              % Nt x 1

        BW = bsxfun(@gt, E, thrE);             % Nt x Nr
        BW(~activeT,:) = false;

        % ---- Enforce TOP-K per active time row (prevents full vertical band fill) ----
        for t = 1:Nt
            if ~activeT(t), continue; end

            idx = find(BW(t,:));
            if isempty(idx)
                % fallback: pick top-K positive excess bins even if threshold missed
                [vals, ord] = sort(E(t,:), 'descend');
                ord = ord(vals > 0);
                if isempty(ord), continue; end
                take = ord(1:min(Kcap, numel(ord)));
                BW(t,:) = false;
                BW(t,take) = true;
            else
                if numel(idx) > Kcap
                    [~, ord] = sort(E(t,idx), 'descend');
                    keep = idx(ord(1:Kcap));
                    row = false(1,Nr);
                    row(keep) = true;
                    BW(t,:) = row;
                end
            end
        end

        % cleanup specks
        BW = bwareaopen(BW, opts.minArea);

        % diagonal connectivity
        if opts.ridgeLen > 1
            se1 = strel('line', opts.ridgeLen, 45);
            if opts.useBothDirs
                se2 = strel('line', opts.ridgeLen, 135);
                BW  = imclose(BW, se1) | imclose(BW, se2);
            else
                BW  = imclose(BW, se1);
            end
        end

        if opts.dilateRad > 0
            BW = imdilate(BW, strel('disk', opts.dilateRad, 0));
        end

        % choose the "best" BW: non-empty but not too full
        nn = nnz(BW);
        fullness = nn / (Nt*Nr + eps);

        % prefer: some pixels, not too many
        score = nn;
        if fullness > 0.35
            score = score - 1e6;  % reject near-full-band solutions
        end
        if nn < 30
            score = score - 1e5;  % reject near-empty
        end

        if score > bestScore
            bestScore = score;
            BWbest = BW;
        end
    end

    % ---- Write back into mask ----
    maskNtNf(:, idxRadarF) = noiseVal;
    tmp = maskNtNf(:, idxRadarF);
    tmp(BWbest) = radarVal;
    maskNtNf(:, idxRadarF) = tmp;
end
