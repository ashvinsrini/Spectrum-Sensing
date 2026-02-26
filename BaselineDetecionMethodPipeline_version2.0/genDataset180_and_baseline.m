%% ===================== DATASET (180) + BASELINE (CYCLO + RADAR) =====================
% Generates ONE consolidated dataset that mixes:
%   (A) 100x one-tech (NR/LTE/WLAN/RADAR)
%   (B)  50x RADAR + one non-RADAR (NR/LTE/WLAN)
%   (C)  30x two non-RADAR with partial overlap probability pOverlap
% and then runs the Baseline_detection_pipeline-style detector on the SAME xMix.
%
% Outputs:
%   dataset_out_180/
%     gt8/   : ground-truth masks with overlap-classes (8 classes)
%     gt5/   : ground-truth masks mapped to 5 classes (overlap -> LEFT tech)
%     baseline5/ : baseline predicted masks (5 classes) + score maps
%     inputfiles_gt8/ : spec_*.png + mask_*.hdf for DL training (8 classes)
%     inputfiles_gt5/ : spec_*.png + mask_*.hdf for baseline evaluation (5 classes)
%
% Required helper functions (you said you already have them):
%   makeTemplatesAtFsCommon, mySpectrogramDB, classId, maskToRGB_values_local,
%   refineRadarMaskAmp_local, widenRadarPulses_local,
%   detectBandsFromSpectrogram, bandEdgesToF0BW, bandpassDDCdecimate,
%   slidingCpCorrCoeffFast, radarScoreFast
%
% NOTE on labels:
%   8-class GT (uint8 values):
%     Noise=0, NR=63, LTE=127, WLAN=191, RADAR=6, NRLTE=32, NRWLAN=96, LTEWLAN=160
%   5-class GT/Baseline (uint8 values):
%     Noise=0, NR=63, LTE=127, WLAN=191, RADAR=6
%
% Save this file as: genDataset180_and_baseline.m and run.

clc; close all; clear;
rng(1);

%% ---------------- User knobs ----------------
N_oneTech   = 100;
N_radarPlus = 50;
N_overlap2  = 30;
Ntotal      = N_oneTech + N_radarPlus + N_overlap2;

% overlap controls for the two-nonRADAR cases
pOverlap   = 0.75;   % probability of partial overlap within the 30 overlap2 trials
minOvFrac  = 0.15;   % min overlap fraction of min(BW1,BW2)
maxOvFrac  = 0.60;   % max overlap fraction of min(BW1,BW2)  (prevents full overlap)

% radar time stripe widening
radarPulseWiden_us = 200;  % widen each radar pulse by +/-200 us

% Output
outRoot = fullfile(pwd,'dataset_out_180');

% Save options
outSize = [256 256];
useFixedCLim = true;
clim_dB = [-40 60];
cmap = parula(256);

%% ---------------- Label maps ----------------
% 5-class values
val5.Noise = uint8(0);
val5.NR    = uint8(63);
val5.LTE   = uint8(127);
val5.WLAN  = uint8(191);
val5.RADAR = uint8(6);
labelVals5 = uint8([0 63 127 191 6]); % index = (trueOfFreq 0..4)+1

% 8-class values (DL GT)
val8 = val5;
val8.NRLTE   = uint8(32);
val8.NRWLAN  = uint8(96);
val8.LTEWLAN = uint8(160);

% For mask RGB visualization (extend your 5 colors with 3 more)
maskColors8_u8 = uint8(round(255 * [ ...
    0.20 0.90 0.90;   % Noise
    0.10 0.60 1.00;   % NR
    0.70 0.40 0.90;   % LTE
    1.00 0.00 1.00;   % WLAN
    1.00 0.65 0.00;   % RADAR
    0.30 0.85 0.10;   % NRLTE
    0.85 0.25 0.10;   % NRWLAN
    0.20 0.20 0.85;   % LTEWLAN
]));

% Class name orders for saving a metadata file / for pixelLabelDatastore
classNames8 = ["Noise","NR","LTE","WLAN","RADAR","NRLTE","NRWLAN","LTEWLAN"];
pixelLabelID8 = double([val8.Noise val8.NR val8.LTE val8.WLAN val8.RADAR val8.NRLTE val8.NRWLAN val8.LTEWLAN]);

classNames5 = ["Noise","NR","LTE","WLAN","RADAR"];
pixelLabelID5 = double([val5.Noise val5.NR val5.LTE val5.WLAN val5.RADAR]);

%% ---------------- Common signal parameters ----------------
FsCommon  = 300e6;
Ttotal_ms = 40;
Nsig      = round(Ttotal_ms*1e-3*FsCommon);

BW.NR    = 20e6;
BW.LTE   = 20e6;
BW.WLAN  = 20e6;
BW.RADAR = 25e6;

fMin = -150e6;  fMax = +150e6;
guardHz = 8e6;

snr_dB.NR    = 10;
snr_dB.LTE   = 10;
snr_dB.WLAN  = 10;
snr_dB.RADAR = 10;

% Radar: force multi regions
numRadarRegions = 4;
radarRandShift  = true;

% Spectrogram params
stftWin  = 4096;
stftHop  = 4096;
stftNfft = 4096;

% Fast processing / decimation (baseline)
FsProc_des = 30e6;
D = max(1, round(FsCommon/FsProc_des));
FsProc = FsCommon / D;

Fcut20 = 12e6;
firOrd = 192;
h20 = fir1(firOrd, Fcut20/(FsCommon/2));
dec20 = dsp.FIRDecimator(D, h20);

Fcut15 = 10e6;
h15 = fir1(firOrd, Fcut15/(FsCommon/2));
dec15 = dsp.FIRDecimator(D, h15);

n = (0:Nsig-1).';

% CP lag definitions (baseline)
FsWLAN_ref   = 20e6;
NfftWLAN_ref = 64;
FsLTE_ref    = 30.72e6;
NfftLTE_ref  = 2048;
carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = 30;
carrier.NSizeGrid         = 51;
carrier.NSlot             = 0;
carrier.NFrame            = 0;
ofdmInfo = nrOFDMInfo(carrier);
FsNR_ref   = ofdmInfo.SampleRate;
NfftNR_ref = ofdmInfo.Nfft;

L_NR_p   = max(1, round((NfftNR_ref  /FsNR_ref  ) * FsProc));
L_LTE_p  = max(1, round((NfftLTE_ref /FsLTE_ref ) * FsProc));
L_WLAN_p = max(1, round((NfftWLAN_ref/FsWLAN_ref) * FsProc));

%% ---------------- Baseline detector knobs ----------------

% Stripe scoring knobs
tauPix   = 0.0;   % pixel threshold used in Baseline_detection_pipeline

%% ---------------- Output folders ----------------
% GT8
gt8.specDir    = fullfile(outRoot,'gt8','spectrogram');     ensureDir_local(gt8.specDir);
gt8.maskDir    = fullfile(outRoot,'gt8','mask');            ensureDir_local(gt8.maskDir);
gt8.maskEnhDir = fullfile(outRoot,'gt8','maskenhanced');    ensureDir_local(gt8.maskEnhDir);
gt8.maskHdfDir = fullfile(outRoot,'gt8','maskhdf');         ensureDir_local(gt8.maskHdfDir);

gt8.inputDir   = fullfile(outRoot,'inputfiles_gt8');        ensureDir_local(gt8.inputDir);

% GT5 (overlap mapped to LEFT)
gt5.specDir    = fullfile(outRoot,'gt5','spectrogram');     ensureDir_local(gt5.specDir);
gt5.maskDir    = fullfile(outRoot,'gt5','mask');            ensureDir_local(gt5.maskDir);
gt5.maskEnhDir = fullfile(outRoot,'gt5','maskenhanced');    ensureDir_local(gt5.maskEnhDir);
gt5.maskHdfDir = fullfile(outRoot,'gt5','maskhdf');         ensureDir_local(gt5.maskHdfDir);

gt5.inputDir   = fullfile(outRoot,'inputfiles_gt5');        ensureDir_local(gt5.inputDir);

% Baseline predicted (5-class)
bl.specDir     = fullfile(outRoot,'baseline5','spectrogram');    ensureDir_local(bl.specDir);
bl.maskDir     = fullfile(outRoot,'baseline5','pred_mask');      ensureDir_local(bl.maskDir);
bl.maskEnhDir  = fullfile(outRoot,'baseline5','pred_maskenh');   ensureDir_local(bl.maskEnhDir);
bl.maskHdfDir  = fullfile(outRoot,'baseline5','pred_maskhdf');   ensureDir_local(bl.maskHdfDir);
bl.matDir      = fullfile(outRoot,'baseline5','mat');            ensureDir_local(bl.matDir);

%% ---------------- Precompute templates once ----------------
tmpl = makeTemplatesAtFsCommon(FsCommon, Nsig, carrier);

%% ---------------- Build mixed schedule (180 trials) ----------------
meta = struct();
meta.Ntotal = Ntotal;
meta.classNames8 = classNames8;
meta.pixelLabelID8 = pixelLabelID8;
meta.classNames5 = classNames5;
meta.pixelLabelID5 = pixelLabelID5;

trialType = strings(Ntotal,1);
trialTechA = strings(Ntotal,1);
trialTechB = strings(Ntotal,1);
trialDoOverlap = false(Ntotal,1);

idx = 1;

% (A) one-tech schedule, balanced
classesA = ["NR","LTE","WLAN","RADAR"];
countsA  = repmat(floor(N_oneTech/numel(classesA)),1,numel(classesA));
remA     = N_oneTech - sum(countsA);
if remA>0
    extra = randperm(numel(classesA), remA);
    countsA(extra) = countsA(extra) + 1;
end
for c=1:numel(classesA)
    k = countsA(c);
    if k>0
        trialType(idx:idx+k-1) = "ONE";
        trialTechA(idx:idx+k-1) = classesA(c);
        idx = idx + k;
    end
end

% (B) radar + one non-radar, balanced
classesB = ["NR","LTE","WLAN"];
countsB  = repmat(floor(N_radarPlus/numel(classesB)),1,numel(classesB));
remB     = N_radarPlus - sum(countsB);
if remB>0
    extra = randperm(numel(classesB), remB);
    countsB(extra) = countsB(extra) + 1;
end
for c=1:numel(classesB)
    k = countsB(c);
    if k>0
        trialType(idx:idx+k-1) = "RADARPLUS";
        trialTechA(idx:idx+k-1) = classesB(c);
        idx = idx + k;
    end
end

% (C) overlap two non-radar, balanced over 3 pairs
pairList = ["NR","LTE";
            "NR","WLAN";
            "LTE","WLAN"];
Npairs = size(pairList,1);
countsC = repmat(floor(N_overlap2/Npairs),1,Npairs);
remC    = N_overlap2 - sum(countsC);
if remC>0
    extra = randperm(Npairs, remC);
    countsC(extra) = countsC(extra) + 1;
end
for p=1:Npairs
    k = countsC(p);
    if k>0
        trialType(idx:idx+k-1)  = "OVERLAP2";
        trialTechA(idx:idx+k-1) = pairList(p,1);
        trialTechB(idx:idx+k-1) = pairList(p,2);
        idx = idx + k;
    end
end

% overlap flags only for OVERLAP2 trials
ovIdx = find(trialType=="OVERLAP2");
Nover = round(pOverlap * numel(ovIdx));
trialDoOverlap(ovIdx(1:Nover)) = true;
trialDoOverlap(ovIdx(Nover+1:end)) = false;

% shuffle all trials together
perm = randperm(Ntotal);
trialType = trialType(perm);
trialTechA = trialTechA(perm);
trialTechB = trialTechB(perm);
trialDoOverlap = trialDoOverlap(perm);

%% ---------------- Main loop ----------------
for i = 1:Ntotal

    typ = trialType(i);
    techA = string(trialTechA(i));
    techB = string(trialTechB(i));
    doOv  = trialDoOverlap(i);

    labelsTrue  = strings(0,1);
    bandsTrueHz = zeros(0,2);
    f0ListHz    = zeros(0,1);

    % ---------------- (A) ONE: one technology ----------------
    if typ == "ONE"
        if techA == "RADAR"
            K = numRadarRegions;
            labelsTrue  = repmat("RADAR", K, 1);
            bandsTrueHz = zeros(K,2);
            f0ListHz    = zeros(K,1);

            BWk   = BW.RADAR;
            f0min = fMin + BWk/2 + guardHz;
            f0max = fMax - BWk/2 - guardHz;

            for r=1:K
                placed=false;
                for attempt=1:300
                    f0 = f0min + (f0max - f0min)*rand();
                    fL = f0 - BWk/2; fH = f0 + BWk/2;
                    if r==1
                        placed=true;
                    else
                        prev = bandsTrueHz(1:r-1,:);
                        if all((fH < prev(:,1)-guardHz) | (fL > prev(:,2)+guardHz))
                            placed=true;
                        end
                    end
                    if placed
                        bandsTrueHz(r,:)=[fL fH];
                        f0ListHz(r)=f0;
                        break;
                    end
                end
                if ~placed
                    error('ONE/RADAR: could not place radar bands.');
                end
            end
        else
            labelsTrue  = techA;
            BWk = BW.(techA);
            f0min = fMin + BWk/2 + guardHz;
            f0max = fMax - BWk/2 - guardHz;
            f0 = f0min + (f0max - f0min)*rand();
            f0ListHz = f0;
            bandsTrueHz = [f0 - BWk/2, f0 + BWk/2];
        end

    % ---------------- (B) RADARPLUS: RADAR regions + one non-radar ----------------
    elseif typ == "RADARPLUS"
        otherTech = techA; % NR/LTE/WLAN

        K = numRadarRegions;
        labelsTrue = [repmat("RADAR",K,1); otherTech];
        bandsTrueHz = zeros(K+1,2);
        f0ListHz    = zeros(K+1,1);

        % Place the other tech first
        BWk = BW.(otherTech);
        f0min = fMin + BWk/2 + guardHz;
        f0max = fMax - BWk/2 - guardHz;
        f0o = f0min + (f0max - f0min)*rand();
        fLo = f0o - BWk/2; fHo = f0o + BWk/2;
        bandsTrueHz(end,:) = [fLo fHo];
        f0ListHz(end)      = f0o;

        % Place radar regions non-overlapping with each other and with other band
        BWkR = BW.RADAR;
        f0minR = fMin + BWkR/2 + guardHz;
        f0maxR = fMax - BWkR/2 - guardHz;

        for r=1:K
            placed=false;
            for attempt=1:500
                f0r = f0minR + (f0maxR - f0minR)*rand();
                fLr = f0r - BWkR/2; fHr = f0r + BWkR/2;

                prev = bandsTrueHz(1:max(0,r-1),:);
                okRadar = true;
                if r>1
                    okRadar = all((fHr < prev(:,1)-guardHz) | (fLr > prev(:,2)+guardHz));
                end
                okOther = (fHr < fLo - guardHz) || (fLr > fHo + guardHz);

                if okRadar && okOther
                    placed=true;
                    bandsTrueHz(r,:) = [fLr fHr];
                    f0ListHz(r)      = f0r;
                    break;
                end
            end
            if ~placed
                error('RADARPLUS: could not place radar region(s) w.r.t. otherTech.');
            end
        end

    % ---------------- (C) OVERLAP2: two non-radar tech with optional partial overlap ----------------
    elseif typ == "OVERLAP2"
        if techB == ""
            error('OVERLAP2 needs techB.');
        end

        labelsTrue  = [techA; techB];
        bandsTrueHz = zeros(2,2);
        f0ListHz    = zeros(2,1);

        % Place techA
        BW1 = BW.(techA);
        f0min1 = fMin + BW1/2 + guardHz;
        f0max1 = fMax - BW1/2 - guardHz;
        f01 = f0min1 + (f0max1 - f0min1)*rand();
        fL1 = f01 - BW1/2; fH1 = f01 + BW1/2;
        bandsTrueHz(1,:) = [fL1 fH1];
        f0ListHz(1) = f01;

        % Place techB relative to techA
        BW2 = BW.(techB);
        f0min2 = fMin + BW2/2 + guardHz;
        f0max2 = fMax - BW2/2 - guardHz;

        placed=false;
        if doOv
            ovMin = max(f0min2, fL1 - BW2/2);
            ovMax = min(f0max2, fH1 + BW2/2);
            minBW = min(BW1,BW2);
            for attempt=1:800
                if ~(ovMin < ovMax), break; end
                f02 = ovMin + (ovMax - ovMin)*rand();
                fL2 = f02 - BW2/2; fH2 = f02 + BW2/2;

                overlapBW = min(fH1,fH2) - max(fL1,fL2);
                if overlapBW <= 0, continue; end

                fullContain12 = (fL2 >= fL1) && (fH2 <= fH1);
                fullContain21 = (fL1 >= fL2) && (fH1 <= fH2);
                if fullContain12 || fullContain21, continue; end

                if overlapBW < minOvFrac*minBW, continue; end
                if overlapBW > maxOvFrac*minBW, continue; end

                placed=true;
                break;
            end
            if ~placed
                doOv=false; % fallback to non-overlap
            end
        end

        if ~doOv
            for attempt=1:800
                f02 = f0min2 + (f0max2 - f0min2)*rand();
                fL2 = f02 - BW2/2; fH2 = f02 + BW2/2;
                if (fH2 < fL1 - guardHz) || (fL2 > fH1 + guardHz)
                    placed=true;
                    break;
                end
            end
        end

        if ~placed
            error('OVERLAP2: could not place techB relative to techA. Reduce guardHz.');
        end

        bandsTrueHz(2,:) = [fL2 fH2];
        f0ListHz(2)      = f02;

    else
        error('Unknown trialType: %s', typ);
    end

    presentNamesUnique = unique(labelsTrue,'stable');

    % ---------------- Build xMix ----------------
    noise = 0.001*(randn(Nsig,1,'single') + 1j*randn(Nsig,1,'single'))/sqrt(2);
    xMix  = noise;

    for kBand = 1:numel(labelsTrue)
        tech = char(labelsTrue(kBand));
        f0   = f0ListHz(kBand);
        x0   = single(tmpl.(tech)(:));

        if strcmpi(tech,'RADAR') && radarRandShift
            sh = randi([0, Nsig-1], 1, 1);
            x0 = circshift(x0, sh);
            x0 = x0 .* exp(1j*2*pi*rand());
        end
        if strcmpi(tech,'RADAR')
            x0 = widenRadarPulses_local(x0, FsCommon, radarPulseWiden_us);
        end

        x0 = x0 ./ (rms(x0)+eps('single')) * sqrt(10^(snr_dB.(tech)/10));
        osc  = exp(1j*2*pi*(f0/FsCommon)*n);
        xMix = xMix + x0 .* single(osc);
    end
    xMix = single(xMix(:));

    % ---------------- Spectrogram ----------------
    [Sdb, fAxis, tAxis] = mySpectrogramDB(xMix, FsCommon, stftWin, stftHop, stftNfft);
    Nf = numel(fAxis);
    Nt = numel(tAxis);

    % ---------------- GT masks (5-class base) ----------------
    % trueOfFreq_base: 0..4 (Noise, NR, LTE, WLAN, RADAR)
    trueOfFreq_base = zeros(Nf,1,'uint8');

    if typ == "OVERLAP2"
        % paint tech1/tech2, and resolve overlap to LEFT tech for 5-class GT
        fL1 = bandsTrueHz(1,1); fH1 = bandsTrueHz(1,2);
        fL2 = bandsTrueHz(2,1); fH2 = bandsTrueHz(2,2);
        tech1 = string(labelsTrue(1));
        tech2 = string(labelsTrue(2));

        idx1 = (fAxis>=fL1) & (fAxis<=fH1);
        idx2 = (fAxis>=fL2) & (fAxis<=fH2);
        idxOv = idx1 & idx2;

        % LEFT by lower start frequency
        if fL1 <= fL2
            leftTech = tech1; rightTech = tech2;
            idxLeftOnly  = idx1 & ~idxOv;
            idxRightOnly = idx2 & ~idxOv;
        else
            leftTech = tech2; rightTech = tech1;
            idxLeftOnly  = idx2 & ~idxOv;
            idxRightOnly = idx1 & ~idxOv;
        end

        idLeft  = uint8(classId(char(leftTech)));
        idRight = uint8(classId(char(rightTech)));

        trueOfFreq_base(idxLeftOnly | idxOv) = idLeft;
        trueOfFreq_base(idxRightOnly)        = idRight;

    else
        % ONE / RADARPLUS: paint each allocated band (no overlap in your placement)
        for kT = 1:size(bandsTrueHz,1)
            fL = bandsTrueHz(kT,1); fH = bandsTrueHz(kT,2);
            id = uint8(classId(char(labelsTrue(kT))));
            idxF = (fAxis >= fL) & (fAxis <= fH);
            trueOfFreq_base(idxF) = id;
        end
    end

    % build Nt x Nf 5-class GT mask values
    maskOfFreq5 = labelVals5(double(trueOfFreq_base)+1);
    gtMask5 = repmat(maskOfFreq5(:).', Nt, 1);

    % If radar exists, time-localize radar pixels (same options you used)
    if any(trueOfFreq_base == uint8(4))
        optsRadarAmp.method         = 'hysteresis';
        optsRadarAmp.highQuantile   = 0.985;
        optsRadarAmp.lowQuantile    = 0.930;
        optsRadarAmp.minArea        = 8;
        optsRadarAmp.diagCloseLen   = 0;
        optsRadarAmp.dilateRad      = 0;
        optsRadarAmp.collapseRuns   = true;
        optsRadarAmp.keepPerRun     = 1;
        optsRadarAmp.maxKeepPerRow  = 4;
        gtMask5 = refineRadarMaskAmp_local(gtMask5, Sdb, trueOfFreq_base, labelVals5, optsRadarAmp);
    end

    % ---------------- GT masks (8-class, overlap pixels become overlap-class) ----------------
    gtMask8 = gtMask5; % start from resolved 5-class mask (already radar-localized)

    if typ == "OVERLAP2"
        % Overwrite only the OVERLAP frequency bins with {NRLTE, NRWLAN, LTEWLAN}
        fL1 = bandsTrueHz(1,1); fH1 = bandsTrueHz(1,2);
        fL2 = bandsTrueHz(2,1); fH2 = bandsTrueHz(2,2);
        tech1 = string(labelsTrue(1));
        tech2 = string(labelsTrue(2));

        idx1 = (fAxis>=fL1) & (fAxis<=fH1);
        idx2 = (fAxis>=fL2) & (fAxis<=fH2);
        idxOv = idx1 & idx2;

        if any(idxOv)
            ovVal = overlapVal8(tech1, tech2, val8);
            gtMask8(:, idxOv) = ovVal;  % constant over time
        end

        % Also ensure single-tech-only regions keep their single-tech values (already true)
    end

    % ---------------- Baseline detection (from Baseline_detection_pipeline) ----------------

    classesBL = ["NR","LTE","WLAN","RADAR"]; 
    Cbl = numel(classesBL);  % 4

    % ===== Baseline_detection_pipeline-style detector on current xMix =====
    % PRI/PW ranges + envelope rate (same as your pipeline)
    PRI_range_s = [1e-3 2e-3];
    PW_range_s  = [50e-6 100e-6];
    FsEnv       = 2e6;

    % (Nwin/hop control the sliding CP-metric smoothing)
    numSymWin = 1; 
    hopFactor = 0.25; 
    symLenUse = 4096;
    Nwin = max(256, round(numSymWin * symLenUse));
    hop  = max(1,   round(hopFactor * symLenUse));

    % 1) Detect frequency bands from spectrogram
    bandsDetHz = detectBandsFromSpectrogram(Sdb, fAxis);
    Kdet = size(bandsDetHz,1);

    stripeOfFreq = zeros(Nf,1,'uint8');
    for b = 1:Kdet
        idxF = (fAxis >= bandsDetHz(b,1)) & (fAxis <= bandsDetHz(b,2));
        stripeOfFreq(idxF) = uint8(b);
    end

    % 2) Stripe scores (NR/LTE/WLAN cyclo + RADAR score)
    stripeScores = zeros(Kdet, Cbl);  % [NR LTE WLAN RADAR]

    for b = 1:Kdet
        fL = bandsDetHz(b,1);
        fH = bandsDetHz(b,2);
        fc = 0.5*(fL+fH);
        bwEst = (fH - fL);

        % choose decimator based on BW estimate (same spirit as your baseline)
        if bwEst <= 16e6
            decObj = dec15;
        else
            decObj = dec20;
        end
        reset(decObj);
        FsBand = FsCommon / D;

        % DDC + decimate
        osc  = exp(-1j*2*pi*(fc/FsCommon)*n);
        xBand = decObj(double(xMix).*osc);

        rhoNR   = slidingCpCorrCoeffFast(xBand, FsBand, L_NR_p,   Nwin, hop);
        rhoLTE  = slidingCpCorrCoeffFast(xBand, FsBand, L_LTE_p,  Nwin, hop);
        rhoWLAN = slidingCpCorrCoeffFast(xBand, FsBand, L_WLAN_p, Nwin, hop);

        sRadar  = radarScoreFast(xBand, FsBand, PRI_range_s, PW_range_s, FsEnv);

        stripeScores(b,:) = [mean(rhoNR), mean(rhoLTE), mean(rhoWLAN), sRadar];
    end

    % 3) Score maps for ROC: Nt x Nf x 4 (constant over time per stripe)
    scoreNtNf = zeros(Nt, Nf, Cbl, 'single');
    bIdx = double(stripeOfFreq(:));
    for c = 1:Cbl
        sFreq = zeros(Nf,1,'single');
        ii = (bIdx > 0);
        if any(ii)
            sFreq(ii) = single(stripeScores(bIdx(ii), c));
        end
        scoreNtNf(:,:,c) = repmat(sFreq(:).', Nt, 1);
    end

    % match orientation of saved images
    scoreNtNf = flipud(scoreNtNf);

    % resize to outSize to align with saved masks
    if ~isempty(outSize)
        scoreMaps = zeros(outSize(1), outSize(2), Cbl, 'single');
        for c = 1:Cbl
            scoreMaps(:,:,c) = imresize(scoreNtNf(:,:,c), outSize, 'bilinear');
        end
    else
        scoreMaps = scoreNtNf;
    end

    % 4) Predicted mask (5-class): argmax stripe + threshold tauPix
    predStripeId = zeros(Kdet,1,'uint8');
    bestStripeScore = zeros(Kdet,1);
    for b = 1:Kdet
        [bestStripeScore(b), idxBest] = max(stripeScores(b,:));
        if bestStripeScore(b) > tauPix
            predStripeId(b) = uint8(idxBest); % 1..4
        else
            predStripeId(b) = uint8(0);
        end
    end

    predOfFreq = zeros(Nf,1,'uint8');
    for jf = 1:Nf
        sb = stripeOfFreq(jf);
        if sb > 0
            predOfFreq(jf) = predStripeId(sb);
        end
    end

    predMask5 = labelVals5(double(predOfFreq)+1);
    predMask5 = repmat(predMask5(:).', Nt, 1);

    % Save ROC inputs for this trial
    save(fullfile(bl.matDir, sprintf('roc_%04d.mat', i)), ...
        'scoreMaps','labelVals5','classesBL','bandsDetHz','stripeScores','stripeOfFreq', ...
        'predStripeId','bestStripeScore','typ','techA','techB','doOv','labelsTrue','bandsTrueHz', ...
        '-v7.3');

% ---------------- Save (GT8, GT5, Baseline5) ----------------
    specRGB = SdbToRGB_local(Sdb, cmap, useFixedCLim, clim_dB);

    % flip to match axis xy
    specRGB_flip = flipud(specRGB);
    %specRGB_flip = specRGB;
    gtMask8_flip = flipud(gtMask8);
    gtMask5_flip = flipud(gtMask5);
    predMask5_flip = flipud(predMask5);

    % Resize
    if ~isempty(outSize)
        specRGB_flip  = imresize(specRGB_flip,  outSize, 'bilinear');
        gtMask8_flip  = imresize(uint8(gtMask8_flip),  outSize, 'nearest');
        gtMask5_flip  = imresize(uint8(gtMask5_flip),  outSize, 'nearest');
        predMask5_flip= imresize(uint8(predMask5_flip),outSize, 'nearest');
    end

    % filenames
    fnSpec = sprintf('spec_%04d.png', i);
    fnM8   = sprintf('mask_%04d.hdf', i);
    fnM5   = sprintf('mask_%04d.hdf', i);

    % --- GT8 ---
    imwrite(specRGB_flip, fullfile(gt8.specDir, fnSpec));
    imwrite(uint8(gtMask8_flip), fullfile(gt8.maskDir, strrep(fnSpec,'spec_','mask_')));

    hdfPath8 = fullfile(gt8.maskHdfDir, sprintf('mask_%04d.hdf', i));
    if exist(hdfPath8,'file'), delete(hdfPath8); end
    imwrite(uint8(gtMask8_flip), hdfPath8, 'hdf');

    maskRGB8 = maskToRGB_values_local(uint8(gtMask8_flip), uint8(pixelLabelID8), maskColors8_u8);
    imwrite(maskRGB8, fullfile(gt8.maskEnhDir, sprintf('maskenh_%04d.png', i)));

    % copy into inputfiles_gt8 (spec png + mask hdf)
    copyfile(fullfile(gt8.specDir, fnSpec), fullfile(gt8.inputDir, fnSpec));
    copyfile(hdfPath8, fullfile(gt8.inputDir, sprintf('mask_%04d.hdf', i)));

    % --- GT5 ---
    imwrite(specRGB_flip, fullfile(gt5.specDir, fnSpec));
    imwrite(uint8(gtMask5_flip), fullfile(gt5.maskDir, strrep(fnSpec,'spec_','mask_')));

    hdfPath5 = fullfile(gt5.maskHdfDir, sprintf('mask_%04d.hdf', i));
    if exist(hdfPath5,'file'), delete(hdfPath5); end
    imwrite(uint8(gtMask5_flip), hdfPath5, 'hdf');

    maskRGB5 = maskToRGB_values_local(uint8(gtMask5_flip), uint8(pixelLabelID5), maskColors8_u8(1:5,:));
    imwrite(maskRGB5, fullfile(gt5.maskEnhDir, sprintf('maskenh_%04d.png', i)));

    copyfile(fullfile(gt5.specDir, fnSpec), fullfile(gt5.inputDir, fnSpec));
    copyfile(hdfPath5, fullfile(gt5.inputDir, sprintf('mask_%04d.hdf', i)));

    % --- Baseline predicted (5-class) ---
    imwrite(specRGB_flip, fullfile(bl.specDir, fnSpec));
    imwrite(uint8(predMask5_flip), fullfile(bl.maskDir, sprintf('predmask_%04d.png', i)));

    hdfPred = fullfile(bl.maskHdfDir, sprintf('predmask_%04d.hdf', i));
    if exist(hdfPred,'file'), delete(hdfPred); end
    imwrite(uint8(predMask5_flip), hdfPred, 'hdf');

    maskPredRGB5 = maskToRGB_values_local(uint8(predMask5_flip), uint8(pixelLabelID5), maskColors8_u8(1:5,:));
    imwrite(maskPredRGB5, fullfile(bl.maskEnhDir, sprintf('predmaskenh_%04d.png', i)));

    save(fullfile(bl.matDir, sprintf('baseline_%04d.mat', i)), ...
        'stripeScores','bandsDetHz','predStripeId','bestStripeScore','stripeOfFreq','typ','techA','techB','doOv', ...
        'bandsTrueHz','labelsTrue');

    fprintf('Saved %3d/%3d | %s | TRUE={%s} | overlap=%d | %s\n', ...
        i, Ntotal, typ, strjoin(presentNamesUnique,','), doOv, fnSpec);

    % store trial metadata in arrays (for later)
    meta.trial(i).type = char(typ);
    meta.trial(i).techA = char(techA);
    meta.trial(i).techB = char(techB);
    meta.trial(i).doOverlap = doOv;
    meta.trial(i).bandsTrueHz = bandsTrueHz;
    meta.trial(i).labelsTrue = labelsTrue;
end

save(fullfile(outRoot,'dataset_meta.mat'), 'meta');

fprintf('\nDONE. Root: %s\n', outRoot);

%% ===================== Local helpers (no extra files needed) =====================
function ensureDir_local(d)
    if ~exist(d,'dir'), mkdir(d); end
end

function RGB = SdbToRGB_local(Sdb, cmap, useFixedCLim, clim_dB)

    % ---- NEW: swap axes for saved image (freq -> x, time -> y) ----
    % mySpectrogramDB usually returns Sdb as [Nf x Nt], so transpose to [Nt x Nf]
    Sdb = Sdb.';   % <--- ADD THIS LINE

    if useFixedCLim
        lo = clim_dB(1); hi = clim_dB(2);
    else
        lo = min(Sdb(:)); hi = max(Sdb(:));
    end
    X = (Sdb - lo) / max(hi-lo, eps);
    X = max(0, min(1, X));
    idx = 1 + floor(X * (size(cmap,1)-1));
    RGB = ind2rgb(idx, cmap);
    RGB = uint8(round(255*RGB));
end


function maskRGB = maskToRGB_values_local(maskU8, labelVals_u8, colors_u8)
% Map a uint8 label mask to RGB given explicit label values.
% labelVals_u8: 1xK uint8 list of label values (e.g., [0 63 127 191 6 ...])
% colors_u8:    Kx3 uint8 RGB colors (same order as labelVals_u8)
    if ~isa(maskU8,'uint8')
        maskU8 = uint8(maskU8);
    end
    labelVals_u8 = uint8(labelVals_u8(:));
    colors_u8 = uint8(colors_u8);

    [H,W] = size(maskU8);
    maskRGB = zeros(H,W,3,'uint8');

    for k = 1:numel(labelVals_u8)
        m = (maskU8 == labelVals_u8(k));
        if any(m(:))
            for ch = 1:3
                tmp = maskRGB(:,:,ch);
                tmp(m) = colors_u8(k,ch);
                maskRGB(:,:,ch) = tmp;
            end
        end
    end
end

function v = overlapVal8(tech1, tech2, val8)
    % Return overlap-class uint8 for the pair (order agnostic)
    a = string(tech1); b = string(tech2);
    if ( (a=="NR"   && b=="LTE")  || (a=="LTE"  && b=="NR") )
        v = val8.NRLTE;
    elseif ( (a=="NR" && b=="WLAN") || (a=="WLAN" && b=="NR") )
        v = val8.NRWLAN;
    elseif ( (a=="LTE" && b=="WLAN") || (a=="WLAN" && b=="LTE") )
        v = val8.LTEWLAN;
    else
        error('Unexpected overlap pair: %s-%s', a, b);
    end
end
