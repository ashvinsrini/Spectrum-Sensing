%% ==============================================================
%  BATCH EXTENSION (10 scenes) ON TOP OF YOUR WORKING SCRIPT
%  - Keeps your per-scene operations the same
%  - Adds ONLY: loop + overlapFlag logic + mask saving + confusion matrix
%
%  Outputs saved to: ./batch_outputs_baseline/
%    - spectrogram_rgb_XX.png
%    - mask_gt_XX.png, mask_gt_rgb_XX.png
%    - mask_baseline_XX.png, mask_baseline_rgb_XX.png
%    - mask_enhanced_XX.png, mask_enhanced_rgb_XX.png   (enhanced = radar time-gated only)
%    - scene_XX.mat (everything stored)
%    - confusion_matrix_rowNorm.png
%% ==============================================================

clear; clc; close all;
rng(7);

%% -------------------- PLOT STYLE (fonts only) --------------------
STYLE.baseFont   = 15;   % tick labels
STYLE.labelFont  = 17;   % xlabel/ylabel
STYLE.titleFont  = 18;   % title
STYLE.legendFont = 14;   % legend
STYLE.cbFont     = 14;   % colorbar ticks
STYLE.textFont   = 16;   % overlay text annotations

set(groot, 'defaultAxesFontWeight','bold', ...
           'defaultTextFontWeight','bold', ...
           'defaultAxesFontSize',STYLE.baseFont, ...
           'defaultTextFontSize',STYLE.baseFont);

%% -------------------- User knobs --------------------
FsCommon   = 100e6;
Ttotal_ms  = 40;

snrNR_dB   = 20;
snrLTE_dB  = 20;
snrWLAN_dB = 20;

% Nominal sub-band centers (Hz)
f0_NR0   = -30e6;
f0_LTE0  =  10e6;
f0_WLAN0 = -30e6;   % NOTE: same as NR -> overlap when both present

BW_NR   = 20e6;
BW_LTE  = 20e6;
BW_WLAN = 20e6;

fade_ms = 0.2;

stftWin  = 4096;
stftHop  = 4096;
stftNfft = 4096;

%% -------------------- Channel fading knobs (already in your code) --------------------
useFading  = true;
DopplerMin = 0; DopplerMax = 500;

%% -------------------- RADAR knobs (UNCHANGED) --------------------
snrRADAR_dB    = 10;
f0_RADAR0      = +30e6;
pulseWidth_us  = 80;
chirpBW_MHz    = 15;
PRI_us         = 1500;
pulsesPerBurst = 2;
numBursts      = 12;

%% -------------------- BATCH knobs (NEW) --------------------
numScenes   = 5;
overlapFlag = false;                 % default false as you asked
jitterHz    = 1e6;                   % +/-1 MHz jitter around the nominal centers

outDir = fullfile(pwd,'batch_outputs_baseline');
if ~exist(outDir,'dir'), mkdir(outDir); end

% Class IDs used in masks:
% 0=Noise, 1=NR, 2=LTE, 3=WLAN, 4=RADAR
classNames = {'LTE','NR','Noise','RADAR','WLAN'};      % order for confusion matrix plot (like your example)
orderIDs   = [2 1 0 4 3];                               % LTE, NR, Noise, RADAR, WLAN

C = zeros(numel(orderIDs), numel(orderIDs));            % accumulate confusion counts over all scenes

%% -------------------- Run batch --------------------
for sceneIdx = 1:numScenes
    fprintf('\n================= SCENE %d / %d =================\n', sceneIdx, numScenes);

    % NEW: decide which comm technologies are present + centers
    sceneCfg = makeSceneConfig(sceneIdx, overlapFlag, f0_NR0, f0_LTE0, f0_WLAN0, f0_RADAR0, jitterHz);

    % NEW: per-scene Doppler (still within your range)
    DopplerHz = DopplerMin + (DopplerMax-DopplerMin)*rand;

    % -------------------- YOUR ORIGINAL PER-SCENE CODE STARTS --------------------
    % (Only tiny additions: include flags + center values from sceneCfg + GT/pred masks + saving)

    f0_NR   = sceneCfg.f0_NR;
    f0_LTE  = sceneCfg.f0_LTE;
    f0_WLAN = sceneCfg.f0_WLAN;
    f0_RADAR = sceneCfg.f0_RADAR;

    includeNR   = sceneCfg.includeNR;
    includeLTE  = sceneCfg.includeLTE;
    includeWLAN = sceneCfg.includeWLAN;

    % Derived lengths
    Ntotal = round(Ttotal_ms*1e-3 * FsCommon);
    Nfade  = round(fade_ms*1e-3 * FsCommon);

    %% =========================================================
    %  1) Generate baseband waveforms at native sampling rates
    %% =========================================================

    % ----- WLAN -----
    cfgWLAN = wlanNonHTConfig('ChannelBandwidth','CBW20','MCS',4,'PSDULength',1000);
    bitsWiFi = randi([0 1], cfgWLAN.PSDULength*8, 1);
    txWLAN   = wlanWaveformGenerator(bitsWiFi, cfgWLAN);
    FsWLAN   = wlanSampleRate(cfgWLAN);

    % ----- NR OFDM -----
    carrier = nrCarrierConfig;
    carrier.SubcarrierSpacing = 30;
    carrier.NSizeGrid         = 106;
    carrier.NSlot             = 0;
    carrier.NFrame            = 0;
    ofdmInfo = nrOFDMInfo(carrier);
    FsNR     = ofdmInfo.SampleRate;

    Nsc  = carrier.NSizeGrid * 12;
    Nsym = carrier.SymbolsPerSlot;
    M    = 16;
    qamSymbols = qammod(randi([0 M-1], Nsc*Nsym, 1), M, 'UnitAveragePower', true);
    gridNR = reshape(qamSymbols, [Nsc Nsym 1]);
    txNR  = nrOFDMModulate(carrier, gridNR);
    txNR  = txNR(:);

    % ----- LTE -----
    if exist('lteRMCDL','file') == 2 && exist('lteRMCDLTool','file') == 2
        rmccfg = lteRMCDL('R.7');
        rmccfg.NCellID = 0;
        rmccfg.NSubframe = 0;
        [txLTE, ~, infoLTE] = lteRMCDLTool(rmccfg, randi([0 1], rmccfg.PDSCH.TrBlkSizes(1), 1));
        FsLTE = infoLTE.SamplingRate;
    else
        FsLTE = 30.72e6;
        txLTE = simpleLTElikeOFDM(FsLTE);
    end

    %% =========================================================
    %  2) Resample all to FsCommon and tile to FULL duration (40 ms)
    %% =========================================================
    xWLAN = resampleAny(txWLAN, FsWLAN, FsCommon);
    xNR   = resampleAny(txNR,   FsNR,   FsCommon);
    xLTE  = resampleAny(txLTE,  FsLTE,  FsCommon);

    xWLAN = repToLen(xWLAN, Ntotal);
    xNR   = repToLen(xNR,   Ntotal);
    xLTE  = repToLen(xLTE,  Ntotal);

    xWLAN = single(xWLAN);
    xNR   = single(xNR);
    xLTE  = single(xLTE);

    %% -------------------- Apply fading at FsCommon (as in your working code) --------------------
    if useFading
        xNR   = applyCommRayleighEPA_FsCommon(xNR,   FsCommon, DopplerHz, 1000 + sceneIdx*10 + 1);
        xLTE  = applyCommRayleighEPA_FsCommon(xLTE,  FsCommon, DopplerHz, 1000 + sceneIdx*10 + 2);
        xWLAN = applyCommRayleighEPA_FsCommon(xWLAN, FsCommon, DopplerHz, 1000 + sceneIdx*10 + 3);
    end

    %% =========================================================
    %  3) Set relative SNRs
    %% =========================================================
    noise = (randn(Ntotal,1,'single') + 1i*randn(Ntotal,1,'single'))/sqrt(2);

    xNR   = xNR   ./ (rms(xNR)   + eps('single')) * sqrt(10^(snrNR_dB/10));
    xLTE  = xLTE  ./ (rms(xLTE)  + eps('single')) * sqrt(10^(snrLTE_dB/10));
    xWLAN = xWLAN ./ (rms(xWLAN) + eps('single')) * sqrt(10^(snrWLAN_dB/10));

    xNR   = applyTaper(xNR,   Nfade);
    xLTE  = applyTaper(xLTE,  Nfade);
    xWLAN = applyTaper(xWLAN, Nfade);

    %% =========================================================
    %  4) Frequency-shift each technology into its own band
    %% =========================================================
    n = single((0:Ntotal-1).');

    xNR_shift   = xNR   .* exp(1i*2*pi*(f0_NR/FsCommon)*n);
    xLTE_shift  = xLTE  .* exp(1i*2*pi*(f0_LTE/FsCommon)*n);
    xWLAN_shift = xWLAN .* exp(1i*2*pi*(f0_WLAN/FsCommon)*n);

    %% =========================================================
    %  RADAR generation (UNCHANGED)
    %% =========================================================
    PW   = pulseWidth_us * 1e-6;
    BW   = chirpBW_MHz   * 1e6;
    PRI  = PRI_us        * 1e-6;

    PW_samp  = max(1, round(PW  * FsCommon));
    PRI_samp = max(PW_samp+1, round(PRI * FsCommon)); 

    Ttotal = Ntotal / FsCommon;

    burstDur = (pulsesPerBurst-1)*PRI + PW;
    tStartMax = max(0, Ttotal - burstDur - 1/FsCommon);
    burstStarts = linspace(0, tStartMax, numBursts);

    xRADAR_bb = zeros(Ntotal,1,'single');
    active = false(Ntotal,1);

    usePhased = (exist('phased.LinearFMWaveform','class') == 8);
    if usePhased
        wfm = phased.LinearFMWaveform( ...
            'SampleRate', FsCommon, ...
            'PulseWidth', PW, ...
            'PRF', 1/PRI, ...
            'SweepBandwidth', BW, ...
            'SweepDirection', 'Up');
    end

    for b = 1:numBursts
        t0 = burstStarts(b);
        for p = 1:pulsesPerBurst
            idx0 = round((t0 + (p-1)*PRI) * FsCommon) + 1;
            idx1 = idx0 + PW_samp - 1;
            if idx0 < 1 || idx1 > Ntotal, continue; end

            if usePhased
                xPRI = wfm(); xPRI = xPRI(:);
                if numel(xPRI) < PW_samp
                    pulse = xPRI;
                else
                    pulse = xPRI(1:PW_samp);
                end
            else
                tt = (0:PW_samp-1).' / FsCommon;
                k = BW / PW;
                f0 = -BW/2;
                pulse = exp(1j*2*pi*( f0*tt + 0.5*k*tt.^2 ));
            end

            w = single(hann(numel(pulse)));
            pulse = single(pulse) .* w;

            xRADAR_bb(idx0:idx1) = xRADAR_bb(idx0:idx1) + pulse;
            active(idx0:idx1) = true;
        end
    end

    if any(active)
        rmsActive = rms(double(xRADAR_bb(active)));
        xRADAR_bb = xRADAR_bb ./ (rmsActive + eps('single')) * sqrt(10^(snrRADAR_dB/10));
    end

    xRADAR_shift = xRADAR_bb .* exp(1i*2*pi*(f0_RADAR/FsCommon)*n);

    %% =========================================================
    %  5) Composite xMix (same structure, just gated by include flags)
    %% =========================================================
    xMix = noise ...
        + (includeNR  * xNR_shift) ...
        + (includeLTE * xLTE_shift) ...
        + (includeWLAN* xWLAN_shift) ...
        + xRADAR_shift;

    %% =========================================================
    %  6) Spectrogram
    %% =========================================================
    [Sdb, f, t] = mySpectrogramDB(xMix, FsCommon, stftWin, stftHop, stftNfft);

    % ---- NEW: Ground truth mask (pixel-wise) ----
    Nt = numel(t); Nf = numel(f);
    maskGT = zeros(Nt, Nf, 'uint8');

    % Comm tech are present full-time when included
    if includeNR
        idxNR = abs(f - f0_NR) <= BW_NR/2;
        maskGT(:, idxNR) = 1;
    end
    if includeLTE
        idxLTE = abs(f - f0_LTE) <= BW_LTE/2;
        maskGT(:, idxLTE) = 2;
    end
    if includeWLAN
        idxWLAN = abs(f - f0_WLAN) <= BW_WLAN/2;
        % NOTE: If overlapFlag=true and NR+WLAN overlap, this single-label GT will overwrite.
        % We keep it simple (single-label mask) because default overlapFlag=false.
        maskGT(:, idxWLAN) = 3;
    end

    % Radar GT is time-gated using the true 'active' samples mapped to STFT frames
    idxRAD = abs(f - f0_RADAR) <= (BW/2);     % radar BW = chirp BW
    radarActiveTF = activeToTimeBins(active, Ntotal, stftWin, stftHop, Nt);
    maskGT(radarActiveTF, idxRAD) = 4;

    % ---- (Your original plots still work; batch saving uses RGB images below) ----

    %% ==========================================================
    %  Robust occupied-band detection (UNCHANGED)
    %% ==========================================================
    P = mean(10.^(Sdb/10), 2);
    f = f(:);

    df = abs(f(2)-f(1));
    smoothHz   = 2e6;
    smoothBins = max(5, round(smoothHz/df));
    Psm = movmean(P, smoothBins);

    idxNoise = Psm <= prctile(Psm, 40);
    muN  = median(Psm(idxNoise));
    sigN = 1.4826 * mad(Psm(idxNoise), 1);

    kOcc   = 6;
    thrOcc = muN + kOcc*sigN;
    occBins0 = (Psm >= thrOcc);

    minBandHz = 8e6;
    gapHz     = 8e6;

    minBandBins = max(1, round(minBandHz/df));
    gapBins     = max(1, round(gapHz/df));

    occBins = fillShortGapsBool(occBins0, gapBins);
    occBins = removeShortRunsBool(occBins,  minBandBins);

    candBandsHz = binsToBands(f, occBins);

    if isempty(candBandsHz)
        warning('No occupied bands detected. Consider lowering kOcc or minBandHz.');
        bandsHz = zeros(0,2);
        stripeLabels = strings(0,1);
    else
        mergeGapHz = 10e6;
        candBandsHz = mergeCloseBands(candBandsHz, mergeGapHz);

        bandPow = zeros(size(candBandsHz,1),1);
        for kk = 1:size(candBandsHz,1)
            idx = (f >= candBandsHz(kk,1)) & (f <= candBandsHz(kk,2));
            bandPow(kk) = sum(Psm(idx));
        end

        minRelPow = 0.15;
        keep = bandPow >= minRelPow*max(bandPow);
        candBandsHz = candBandsHz(keep,:);
        bandPow     = bandPow(keep);

        Kwant = 3;
        [~,ord] = sort(bandPow, 'descend');
        ord = ord(1:min(Kwant, numel(ord)));
        bandsHz = candBandsHz(ord,:);
        bandsHz = sortrows(bandsHz, 1);
    end

    %% ==========================================================
    %  Stripe classification (UNCHANGED core logic)
    %% ==========================================================
    NfftWLAN = 64;
    L_WLAN   = round(NfftWLAN * (FsCommon/FsWLAN));

    NfftNR = ofdmInfo.Nfft;
    L_NR   = round(NfftNR * (FsCommon/FsNR));

    NfftLTE = 1024;
    L_LTE   = round(NfftLTE * (FsCommon/FsLTE));

    techList = {'NR','LTE','WLAN'};
    lagList  = [L_NR, L_LTE, L_WLAN];

    maxLag = max(lagList);
    Nwin   = max(8*maxLag, round(0.5e-3*FsCommon));
    hop    = round(0.25*Nwin);
    firOrd = 256;

    K = size(bandsHz,1);

    stripeLabels = strings(K,1);
    stripeScores = zeros(K, numel(techList)); 
    stripeRadar  = zeros(K,1);                
    PRI_range_s = [1e-3 2e-3];
    PW_range_s  = [50e-6 100e-6];
    tauRadar    = 0.25;
    dutyMax     = 0.20;

    % NEW: baseline mask from predicted bands/labels
    maskPred = zeros(Nt, Nf, 'uint8');

    % NEW: enhanced mask (radar time-gated for nicer visualization)
    maskEnh = zeros(Nt, Nf, 'uint8');

    radarBandIdx = []; radarBandMask = false(Nt,1);

    for kk = 1:K
        band = bandsHz(kk,:);

        [xBand, fcHz] = extractStripeTimeSeries(xMix, FsCommon, band(1), band(2), firOrd); 

        scores = zeros(1,numel(techList));
        for j = 1:numel(techList)
            L = lagList(j);
            [rhoVec, ~] = slidingCpCorrCoeffFast(xBand, FsCommon, L, Nwin, hop);
            scores(j) = mean(rhoVec);
        end

        % Radar test
        FsEnv = 2e6;
        D = max(1, round(FsCommon/FsEnv));
        FsD = FsCommon / D;

        e = abs(xBand).^2;
        eSm = movmean(e, D);
        eD  = eSm(1:D:end);
        xBandD = sqrt(eD);

        [scoreRadar, ~, duty, ~] = radarPulseTrainScore(xBandD, FsD, PRI_range_s, PW_range_s);
        stripeRadar(kk) = scoreRadar;

        if (scoreRadar >= tauRadar) && (duty <= dutyMax)
            stripeLabels(kk) = "RADAR";
        else
            [~,idxMax] = max(scores);
            stripeLabels(kk) = techList{idxMax};
        end

        % -------- NEW: fill predicted masks using detected bands --------
        idxF = (f >= band(1)) & (f <= band(2));
        cls = labelToID(stripeLabels(kk)); % 0..4
        maskPred(:, idxF) = cls;

        % enhanced starts as baseline
        maskEnh(:, idxF) = cls;

        % if radar band, estimate time activity from the extracted stripe (receiver-side)
        if stripeLabels(kk) == "RADAR"
            radarBandIdx = idxF;
            radarBandMask = radarActivityFromBand(xBand, Ntotal, stftWin, stftHop, Nt);
        end
    end

    % NEW: apply radar time-gating ONLY in enhanced mask
    if ~isempty(radarBandIdx)
        maskEnh(~radarBandMask, radarBandIdx) = 0;
    end

    %% ===================== NEW: Save spectrogram + masks (RGB) =====================
    specRGB = spectrogramToRGB(Sdb, true);                % true => flip to match axis xy visual
    gtRGB   = maskToRGB(maskGT, true);
    prRGB   = maskToRGB(maskPred, true);
    enhRGB  = maskToRGB(maskEnh, true);

    imwrite(specRGB, fullfile(outDir, sprintf('spectrogram_rgb_%02d.png', sceneIdx)));
    imwrite(uint8(flipud(maskGT)),   fullfile(outDir, sprintf('mask_gt_%02d.png', sceneIdx)));
    imwrite(gtRGB,  fullfile(outDir, sprintf('mask_gt_rgb_%02d.png', sceneIdx)));

    imwrite(uint8(flipud(maskPred)), fullfile(outDir, sprintf('mask_baseline_%02d.png', sceneIdx)));
    imwrite(prRGB,  fullfile(outDir, sprintf('mask_baseline_rgb_%02d.png', sceneIdx)));

    imwrite(uint8(flipud(maskEnh)),  fullfile(outDir, sprintf('mask_enhanced_%02d.png', sceneIdx)));
    imwrite(enhRGB, fullfile(outDir, sprintf('mask_enhanced_rgb_%02d.png', sceneIdx)));

    save(fullfile(outDir, sprintf('scene_%02d.mat', sceneIdx)), ...
        'sceneCfg','Sdb','f','t','maskGT','maskPred','maskEnh','specRGB', ...
        'snrNR_dB','snrLTE_dB','snrWLAN_dB','snrRADAR_dB','DopplerHz');

  
    % ===================== Robust gtIdx + scoreMaps =====================
    iScene = sceneIdx;   % use your loop index
    
    % --- Decide GT label convention (0..K-1 or 1..K) ---
    minL = double(min(maskGT(:)));
    maxL = double(max(maskGT(:)));
    
    if minL == 0
        gtIdx = uint16(maskGT(:)) + 1;   % -> 1..(maxL+1)
        K = maxL + 1;
    elseif minL == 1
        gtIdx = uint16(maskGT(:));       % already 1..maxL
        K = maxL;
    else
        error('Unexpected maskGT labels: min=%g max=%g', minL, maxL);
    end
    
    % --- Define classNames to match K ---
    if K == 4
        classNames = {'Noise','NR','LTE','WLAN'};
    elseif K == 5
        classNames = {'Noise','NR','LTE','WLAN','RADAR'};
    else
        error('Unsupported number of classes K=%d (maskGT max=%g, min=%g).', K, maxL, minL);
    end
    
    % --- Build score maps with matching number of classes ---
    scoreMaps = buildScoreMapsFromSpectrogramAndMask_K(Sdb, maskEnh, K);
    
    % --- Final sanity check before save ---
    if any(gtIdx < 1) || any(gtIdx > K)
        u = unique(gtIdx);
        error('gtIdx values out of range. K=%d, unique(gtIdx)=%s', K, mat2str(u(:).'));
    end
    
    save(fullfile(outDir, sprintf('scene_%02d_scores.mat', iScene)), ...
         'gtIdx','scoreMaps','classNames','-v7.3');
    % ============================================================================%
        

    %% ===================== NEW: Accumulate confusion (baseline mask only) =====================
    C = C + confusionFromMasks(maskGT, maskPred, orderIDs);

    close all;
end

%%============================ ROC Curves====================================

%% ==========================================================
%  Pixel-level ROC (baseline) aggregated over all saved scenes
%  Requires per-pixel scoreMaps for ROC (H x W x K)
%% ==========================================================

% --- Where your per-scene .mat files are ---
% outDir must already exist in your main script
scoreFiles = dir(fullfile(outDir, 'scene_*_scores.mat'));

if isempty(scoreFiles)
    error(['No scene_*_scores.mat found in outDir. ', ...
           'You need to save gtIdx + scoreMaps per scene (see save line above).']);
end

% --- Accumulate pixels from all scenes (per class) ---
allScoresByClass = [];   % (Npix_total x K)
allGtByClass     = [];   % (Npix_total x 1) class index [1..K] or [0..K-1], we fix below
classNamesGlobal = [];

for iF = 1:numel(scoreFiles)
    D = load(fullfile(scoreFiles(iF).folder, scoreFiles(iF).name));

    if ~isfield(D,'gtIdx') || ~isfield(D,'scoreMaps')
        error('File %s is missing gtIdx or scoreMaps.', scoreFiles(iF).name);
    end

    gtIdx = D.gtIdx;
    scoreMaps = D.scoreMaps;

    if isempty(classNamesGlobal)
        if isfield(D,'classNames')
            classNamesGlobal = D.classNames;
        else
            % fallback names
            K = size(scoreMaps,3);
            classNamesGlobal = "C" + string(1:K);
        end
    end

    % ---- Force gt to be 1..K indexing ----
    K = size(scoreMaps,3);

    gt = double(gtIdx);
    if min(gt(:)) == 0
        gt = gt + 1;  % convert 0..K-1 -> 1..K
    end

    % basic sanity
    if any(gt(:) < 1) || any(gt(:) > K)
        error('gtIdx values out of range in %s. Expected 0..K-1 or 1..K.', scoreFiles(iF).name);
    end

    % flatten pixels
    S = reshape(scoreMaps, [], K);
    G = reshape(gt, [], 1);

    allScoresByClass = [allScoresByClass; S];
    allGtByClass     = [allGtByClass; G]; 
end

classNamesGlobal = string(classNamesGlobal(:).');  % row string array
K = numel(classNamesGlobal);

% --- Compute per-class ROC (one-vs-rest) ---
roc = struct();
roc.fpr = cell(K,1);
roc.tpr = cell(K,1);
roc.auc = zeros(K,1);

for k = 1:K
    yTrue = (allGtByClass == k);       % positives are class k
    yScore = allScoresByClass(:,k);    % score for class k

    [fpr,tpr,auc] = localPerfcurve(yTrue, yScore);
    roc.fpr{k} = fpr;
    roc.tpr{k} = tpr;
    roc.auc(k) = auc;
end

% --- Micro-average ROC (stack all one-vs-rest decisions) ---
yTrueMicro  = false(numel(allGtByClass)*K,1);
yScoreMicro = zeros(numel(allGtByClass)*K,1);

idx0 = 0;
for k = 1:K
    idx1 = idx0 + numel(allGtByClass);
    yTrueMicro(idx0+1:idx1)  = (allGtByClass == k);
    yScoreMicro(idx0+1:idx1) = allScoresByClass(:,k);
    idx0 = idx1;
end
[fprMicro,tprMicro,aucMicro] = localPerfcurve(yTrueMicro, yScoreMicro);

% --- Macro-average ROC (average TPR on common FPR grid) ---
fprGrid = linspace(0,1,501).';
tprInterp = zeros(numel(fprGrid), K);
for k = 1:K
    %tprInterp(:,k) = interp1(roc.fpr{k}, roc.tpr{k}, fprGrid, 'linear', 'extrap');
    [fprU, tprU] = makeUniqueRocXY(roc.fpr{k}, roc.tpr{k});

    if numel(fprU) < 2
        % degenerate ROC curve (all scores same etc.)
        tprInterp(:,k) = tprU(1) * ones(size(fprGrid));
    else
        tprInterp(:,k) = interp1(fprU, tprU, fprGrid, 'linear', 'extrap');
    end
    tprInterp(:,k) = max(0, min(1, tprInterp(:,k)));
end
tprMacro = mean(tprInterp, 2);
aucMacro = trapz(fprGrid, tprMacro);

% --- Plot ROC figure ---
figure('Color','w'); hold on; grid on;
for k = 1:K
    plot(roc.fpr{k}, roc.tpr{k}, 'LineWidth', 1.4);
end
plot(fprMicro, tprMicro, 'k--', 'LineWidth', 2.0);
plot(fprGrid, tprMacro, 'k:',  'LineWidth', 2.2);

xlabel('False positive rate','FontWeight','bold');
ylabel('True positive rate','FontWeight','bold');
title('Pixel-level ROC (baseline scores) aggregated over scenes','FontWeight','bold');

leg = strings(1,K+2);
for k = 1:K
    leg(k) = sprintf('%s (AUC=%.3f)', classNamesGlobal(k), roc.auc(k));
end
leg(K+1) = sprintf('Micro (AUC=%.3f)', aucMicro);
leg(K+2) = sprintf('Macro (AUC=%.3f)', aucMacro);
legend(leg, 'Location','southeast');

set(gca,'FontWeight','bold');

saveas(gcf, fullfile(outDir,'roc_baseline_pixelLevel.png'));

fprintf('\nSaved: %s\n', fullfile(outDir,'roc_baseline_pixelLevel.png'));

%% ===================== NEW helper functions =====================
function [fpr, tpr, auc] = localPerfcurve(yTrue, yScore)
% yTrue: logical column vector
% yScore: numeric column vector (higher => more positive)

    yTrue  = yTrue(:) > 0;
    yScore = double(yScore(:));

    % Use perfcurve if available, else fallback
    if exist('perfcurve','file') == 2
        [fpr,tpr,~,auc] = perfcurve(yTrue, yScore, true);
        fpr = fpr(:); tpr = tpr(:);
        return;
    end

    % Fallback ROC (no toolbox)
    [s,ord] = sort(yScore,'descend');
    y = yTrue(ord);

    P = sum(y==1);
    N = sum(y==0);
    if P==0 || N==0
        fpr = [0;1]; tpr = [0;1];
        auc = NaN;
        return;
    end

    tp = cumsum(y==1);
    fp = cumsum(y==0);

    tpr = tp / P;
    fpr = fp / N;

    % add endpoints
    fpr = [0; fpr; 1];
    tpr = [0; tpr; 1];

    auc = trapz(fpr, tpr);
end

%%============================================================




%% ===================== NEW: Plot/save row-normalized confusion =====================
Cn = rowNormalize(C);

figure('Color','w');
imagesc(100*Cn); axis image;
colormap(parula); colorbar;
title('Confusion matrix (Baseline mask) - row normalized', ...
    'FontSize',STYLE.titleFont,'FontWeight','bold');

set(gca,'XTick',1:numel(classNames),'YTick',1:numel(classNames), ...
    'XTickLabel',classNames,'YTickLabel',classNames, ...
    'FontSize',STYLE.baseFont,'FontWeight','bold', ...
    'XTickLabelRotation',35);

xlabel('Predicted','FontSize',STYLE.labelFont,'FontWeight','bold');
ylabel('True','FontSize',STYLE.labelFont,'FontWeight','bold');

% annotate %
for i = 1:size(Cn,1)
    for j = 1:size(Cn,2)
        txt = sprintf('%.1f%%', 100*Cn(i,j));
        text(j,i,txt,'HorizontalAlignment','center','FontWeight','bold','Color','k');
    end
end

saveas(gcf, fullfile(outDir,'confusion_matrix_rowNorm.png'));

fprintf('\nDone. Outputs saved in:\n  %s\n', outDir);

%% ===================== NEW helper functions (ONLY ADDITIONS) =====================

function sceneCfg = makeSceneConfig(sceneIdx, overlapFlag, f0_NR0, f0_LTE0, f0_WLAN0, f0_RADAR0, jitterHz)
% overlapFlag=false (default):
%   - keep 3 stripes total (for your Kwant=3): LTE + (NR or WLAN) + RADAR
% overlapFlag=true:
%   - LTE + NR + WLAN (NR and WLAN overlap due to same nominal center) + RADAR -> still ~3 bands

    jit = @( ) (2*rand-1)*jitterHz;

    sceneCfg.includeLTE = true;
    sceneCfg.includeNR  = true;
    sceneCfg.includeWLAN= true;

    if ~overlapFlag
        % Alternate between {LTE+NR} and {LTE+WLAN} to avoid NR/WLAN overlap
        if mod(sceneIdx,2)==1
            sceneCfg.includeNR   = true;
            sceneCfg.includeWLAN = false;
        else
            sceneCfg.includeNR   = false;
            sceneCfg.includeWLAN = true;
        end
    end

    sceneCfg.f0_NR    = f0_NR0   + jit();
    sceneCfg.f0_LTE   = f0_LTE0  + jit();
    sceneCfg.f0_WLAN  = f0_WLAN0 + jit();
    sceneCfg.f0_RADAR = f0_RADAR0+ jit();
end

function id = labelToID(lbl)
% Maps stripe label -> mask class ID
    switch string(lbl)
        case "NR",    id = uint8(1);
        case "LTE",   id = uint8(2);
        case "WLAN",  id = uint8(3);
        case "RADAR", id = uint8(4);
        otherwise,    id = uint8(0);
    end
end

function tf = activeToTimeBins(active, Ntotal, win, hop, Nt)
% Maps sample-level active[] -> STFT time-bin active (any active in the window)
    tf = false(Nt,1);
    for k = 1:Nt
        i0 = (k-1)*hop + 1;
        i1 = min(Ntotal, i0 + win - 1);
        tf(k) = any(active(i0:i1));
    end
end

function tf = radarActivityFromBand(xBand, Ntotal, win, hop, Nt)
% Receiver-side radar activity estimate from band-isolated time series.
% Uses window energy + robust threshold (median + 6*MAD).
    E = zeros(Nt,1);
    for k = 1:Nt
        i0 = (k-1)*hop + 1;
        i1 = min(Ntotal, i0 + win - 1);
        seg = xBand(i0:i1);
        E(k) = mean(abs(seg).^2);
    end
    mu = median(E);
    s  = 1.4826*mad(E,1);
    thr = mu + 6*s;
    tf = (E >= thr);
end

function rgb = spectrogramToRGB(Sdb, doFlip)
% Converts Sdb (Nf x Nt) to an RGB image matching imagesc(f,t,Sdb.') with axis xy.
    A = Sdb.';                       % Nt x Nf
    if doFlip, A = flipud(A); end

    lo = prctile(A(:), 5);
    hi = prctile(A(:), 99);
    A = (A - lo) / max(hi-lo, 1e-12);
    A = min(max(A,0),1);

    if exist('turbo','file') == 2
        cmap = turbo(256);
    else
        cmap = parula(256);
    end

    idx = uint8(round(A*255) + 1);   % 1..256
    rgb = ind2rgb(idx, cmap);
    rgb = uint8(round(255*rgb));
end

function rgb = maskToRGB(mask, doFlip)
% Class colors for visibility (reference):
% 0 Noise=black, 1 NR=red, 2 LTE=green, 3 WLAN=yellow, 4 RADAR=cyan
    if doFlip, mask = flipud(mask); end
    cmap = [0 0 0;
            1 0 0;
            0 1 0;
            1 1 0;
            0 1 1];
    mask = uint8(mask);
    mask(mask>4) = 0;
    rgb = ind2rgb(mask+1, cmap);
    rgb = uint8(round(255*rgb));
end

function C = confusionFromMasks(gt, pr, orderIDs)
% Accumulate confusion counts for classes in orderIDs (e.g., [2 1 0 4 3])
    gt = uint8(gt(:));
    pr = uint8(pr(:));
    K = numel(orderIDs);
    C = zeros(K,K);

    for i = 1:K
        gi = orderIDs(i);
        idx = (gt == gi);
        if ~any(idx), continue; end
        pr_i = pr(idx);
        for j = 1:K
            pj = orderIDs(j);
            C(i,j) = C(i,j) + sum(pr_i == pj);
        end
    end
end

function Cn = rowNormalize(C)
% Row-normalize confusion matrix
    C = double(C);
    rs = sum(C,2);
    rs(rs==0) = 1;
    Cn = C ./ rs;
end


function y = applyCommRayleighEPA_FsCommon(x, Fs, dopplerHz, seed)
% Applies a lightweight multipath Rayleigh fading (EPA taps) at FsCommon.
% You still do SNR normalization later, so this is safe for your pipeline.

    x = x(:);

    % Use comm.RayleighChannel if available
    if exist('comm.RayleighChannel','class') == 8
        [pd, pg] = epaTaps_local();

        ch = comm.RayleighChannel( ...
            'SampleRate', Fs, ...
            'PathDelays', pd, ...
            'AveragePathGains', pg, ...
            'MaximumDopplerShift', dopplerHz, ...
            'NormalizePathGains', true, ...
            'RandomStream', 'mt19937ar with seed', ...
            'Seed', seed);

        y = ch(x);
        y = cast(y, 'like', x);
        return;
    end

    % Fallback if Communications Toolbox channel object is not available
    y = applyFlatBlockRayleighFallback(x, Fs, dopplerHz, seed);
    y = cast(y, 'like', x);
end

function [pd, pg] = epaTaps_local()
% EPA delays/gains (seconds / dB)
    pd = [0 30 70 90 110 190 410]*1e-9;
    pg = [0 -1 -2 -3 -8 -17.2 -20.8];
end

function y = applyFlatBlockRayleighFallback(x, Fs, dopplerHz, seed)
% Simple time-varying flat Rayleigh fading fallback.
% Good enough to test robustness if comm.RayleighChannel is missing.

    rng(seed);
    x = x(:);

    if dopplerHz <= 0
        h = (randn + 1j*randn)/sqrt(2);
        y = x*h;
        return;
    end

    % coherence time ~ 0.423/fD, update every ~0.25 Tc (min 0.1 ms)
    Tc = 0.423 / dopplerHz;
    Tblk = max(0.1e-3, 0.25*Tc);
    L = max(2048, round(Tblk*Fs));

    N = numel(x);
    nb = ceil(N/L);

    a = exp(-Tblk / max(Tc,1e-6));
    h = complex(zeros(nb,1));
    h(1) = (randn+1j*randn)/sqrt(2);
    for k = 2:nb
        w = (randn+1j*randn)/sqrt(2);
        h(k) = a*h(k-1) + sqrt(max(0,1-a^2))*w;
    end

    y = x;
    for k = 1:nb
        i0 = (k-1)*L + 1;
        i1 = min(N, k*L);
        y(i0:i1) = y(i0:i1) * h(k);
    end
end


function scoreMaps = buildScoreMapsFromSpectrogramAndMask_K(Sdb, maskPred, K)
% Sdb: (Nf x Nt) in dB
% maskPred: (Nt x Nf) labels {0,1,2,3,(4)}
% scoreMaps: (Nt x Nf x K) in [0,1], channels aligned with classNames

    Sd = double(Sdb.');  % -> (Nt x Nf)

    lo = prctile(Sd(:), 10);
    hi = prctile(Sd(:), 99);
    conf = (Sd - lo) / max(hi - lo, eps);
    conf = single(min(max(conf, 0), 1));

    [Nt, Nf] = size(maskPred);
    scoreMaps = zeros(Nt, Nf, K, 'single');

    % Scores for tech classes (use conf gated by predicted region)
    if K >= 2, scoreMaps(:,:,2) = conf .* single(maskPred == 1); end % NR
    if K >= 3, scoreMaps(:,:,3) = conf .* single(maskPred == 2); end % LTE
    if K >= 4, scoreMaps(:,:,4) = conf .* single(maskPred == 3); end % WLAN
    if K >= 5, scoreMaps(:,:,5) = conf .* single(maskPred == 4); end % RADAR (if used)

    % Noise score as complement of max non-noise score
    if K >= 2
        sig = max(scoreMaps(:,:,2:K), [], 3);
        scoreMaps(:,:,1) = 1 - sig;
    else
        scoreMaps(:,:,1) = 1;
    end
end



function [fprU, tprU] = makeUniqueRocXY(fpr, tpr)
% Ensures FPR is strictly nondecreasing and UNIQUE for interp1.
% For duplicate FPR values, we keep the BEST TPR (max), which is ROC-consistent.

    fpr = double(fpr(:));
    tpr = double(tpr(:));

    % remove NaN/Inf
    ok = isfinite(fpr) & isfinite(tpr);
    fpr = fpr(ok);
    tpr = tpr(ok);

    % clamp to [0,1]
    fpr = min(max(fpr,0),1);
    tpr = min(max(tpr,0),1);

    % sort by fpr
    [fprS, ord] = sort(fpr, 'ascend');
    tprS = tpr(ord);

    % unique FPR; for duplicates take max TPR (upper envelope)
    [fprU, ~, g] = unique(fprS, 'stable');
    tprU = accumarray(g, tprS, [], @max);

    % ensure endpoints exist (helps stability)
    if isempty(fprU)
        fprU = 0; tprU = 0;
        return;
    end
    if fprU(1) > 0
        fprU = [0; fprU];
        tprU = [0; tprU];
    end
    if fprU(end) < 1
        fprU = [fprU; 1];
        tprU = [tprU; tprU(end)];   % keep last TPR
    end
end