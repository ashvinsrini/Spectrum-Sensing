function maskNtNf = refineRadarMaskAmp_local(maskNtNf, Sdb, trueOfFreq, labelVals, opts)
% Amplitude thresholding inside radar band + SAFE removal of horizontal artifacts.
% Never wipes everything (has fallback).

    if nargin < 5 || isempty(opts), opts = struct; end
    if ~isfield(opts,'method'),        opts.method = 'hysteresis'; end
    if ~isfield(opts,'highQuantile'),  opts.highQuantile = 0.985; end
    if ~isfield(opts,'lowQuantile'),   opts.lowQuantile  = 0.930; end
    if ~isfield(opts,'minArea'),       opts.minArea = 8; end
    if ~isfield(opts,'diagCloseLen'),  opts.diagCloseLen = 7; end
    if ~isfield(opts,'dilateRad'),     opts.dilateRad = 1; end

    % --- horizontal pruning knobs (safe defaults) ---
    if ~isfield(opts,'horizMaxH'),     opts.horizMaxH = 2; end
    if ~isfield(opts,'horizMinW'),     opts.horizMinW = 8; end
    if ~isfield(opts,'horizAspect'),   opts.horizAspect = 5; end
    if ~isfield(opts,'horizMaxAng'),   opts.horizMaxAng = 20; end
    if ~isfield(opts,'pruneKeepFrac'), opts.pruneKeepFrac = 0.25; end

    radarId  = uint8(4);
    noiseVal = uint8(labelVals(1));
    radarVal = uint8(labelVals(5));

    Nf = numel(trueOfFreq);

    % ---- Make spectrogram view Stf as Nt x Nf ----
    if size(Sdb,1) == Nf
        Stf = Sdb.';  % (Nf x Nt) -> (Nt x Nf)
    else
        Stf = Sdb;    % (Nt x Nf)
    end
    [Nt, Nf_check] = size(Stf);
    if Nf_check ~= Nf
        error('Mismatch: trueOfFreq=%d bins, spectrogram=%d bins.', Nf, Nf_check);
    end

    idxRadarF = (trueOfFreq(:) == radarId);
    if ~any(idxRadarF), return; end

    Srad = Stf(:, idxRadarF);     % Nt x Nr

    % ---- Robust baseline removal ----
    idxNonRadarF = ~idxRadarF;
    if any(idxNonRadarF)
        base_t = median(Stf(:, idxNonRadarF), 2);  % Nt x 1
    else
        base_t = median(Srad, 2);
    end
    E = bsxfun(@minus, Srad, base_t);     % excess dB
    E = movmedian(E, 3, 2);               % light smoothing along frequency

    % ---- Thresholds from distribution ----
    v = E(:);
    v = v(isfinite(v));
    if isempty(v), return; end

    Thigh = quantile(v, opts.highQuantile);
    Tlow  = quantile(v, opts.lowQuantile);
    if Tlow > Thigh
        tmp = Tlow; Tlow = Thigh; Thigh = tmp;
    end

    % ---- Binary detection ----
    switch lower(opts.method)
        case 'single'
            BW = (E > Tlow);

        otherwise % hysteresis
            seed = (E > Thigh);
            mask = (E > Tlow);
            try
                BW = imreconstruct(seed, mask);
            catch
                BW = mask;
            end
    end

    % ---- Cleanup ----
    BW = bwareaopen(BW, opts.minArea);

    % Optional diagonal closing to connect dotted chirps a bit BEFORE pruning
    if opts.diagCloseLen > 1
        se1 = strel('line', opts.diagCloseLen, 45);
        se2 = strel('line', opts.diagCloseLen, 135);
        BW  = imclose(BW, se1) | imclose(BW, se2);
    end

    BW0 = BW;  % keep a copy (for safety fallback)

    % ============================================================
    % SAFE horizontal artifact removal (component-wise)
    % Remove only components that look like thin horizontal bars.
    % ============================================================
    CC = bwconncomp(BW);
    if CC.NumObjects > 0
        stats = regionprops(CC, 'BoundingBox', 'Orientation');
        kill = false(CC.NumObjects,1);

        for k = 1:CC.NumObjects
            bb = stats(k).BoundingBox;  % [x y w h]
            w = bb(3); h = bb(4);
            asp = w / max(h, 1e-6);

            ang = stats(k).Orientation; % 0 ~ horizontal (MATLAB convention)
            if isnan(ang)
                angOk = true; % don't kill just because angle is NaN (common for tiny blobs)
            else
                angOk = (abs(ang) <= opts.horizMaxAng);
            end

            % Horizontal artifact signature: very thin, long, high aspect ratio, near-horizontal
            if (h <= opts.horizMaxH) && (w >= opts.horizMinW) && (asp >= opts.horizAspect) && angOk
                kill(k) = true;
            end
        end

        BW2 = BW;
        for k = find(kill).'
            BW2(CC.PixelIdxList{k}) = false;
        end
        BW = BW2;
    end

    % Safety fallback: if pruning killed too much, revert
    if nnz(BW0) > 0 && (nnz(BW) < opts.pruneKeepFrac * nnz(BW0))
        BW = BW0;
    end

    % Final small dilation (optional)
    if opts.dilateRad > 0
        BW = imdilate(BW, strel('disk', opts.dilateRad, 0));
    end

    % ---- Write back into full mask ----
    maskNtNf(:, idxRadarF) = noiseVal;
    tmp = maskNtNf(:, idxRadarF);
    tmp(BW) = radarVal;
    maskNtNf(:, idxRadarF) = tmp;
end
