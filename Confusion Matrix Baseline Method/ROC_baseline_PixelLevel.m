%% ===================== PIXEL-LEVEL ROC + CONFUSION (FIXED: stripe region => predicted pixels) =====================
clc; close all; clear;
rng(1);

%% -------- Monte Carlo knobs --------
Nmc = 100;
classes = ["NR","LTE","WLAN","RADAR"];
C = numel(classes);

pPresent   = 0.60;
minNumTech = 1;

NsampPerTrial = 200000;     % sampled pixels per trial

%% -------- Wideband capture --------
FsCommon  = 300e6;
Ttotal_ms = 40;
Ntotal    = round(Ttotal_ms*1e-3*FsCommon);

BW.NR    = 20e6;
BW.LTE   = 20e6;
BW.WLAN  = 20e6;
BW.RADAR = 15e6;

fMin = -150e6;  fMax = +150e6;
guardHz = 8e6;

snr_dB.NR    = 10;
snr_dB.LTE   = 10;
snr_dB.WLAN  = 10;
snr_dB.RADAR = 10;

%% -------- Spectrogram params --------
doPrintDetailTrial = false;
doPlotSpec = false;
plotEvery  = 1;

stftWin  = 4096;
stftHop  = 4096;
stftNfft = 4096;

%% -------- FAST per-stripe processing: DDC + decimate --------
FsProc_des = 30e6;
D = max(1, round(FsCommon/FsProc_des));
FsProc = FsCommon / D;
fprintf('FsCommon=%.0f MHz, FsProc=%.2f MHz (D=%d), Ntotal=%d\n', FsCommon/1e6, FsProc/1e6, D, Ntotal);

Fcut20 = 12e6;
Fcut15 = 10e6;
firOrd = 192;

h20 = fir1(firOrd, Fcut20/(FsCommon/2));
h15 = fir1(firOrd, Fcut15/(FsCommon/2));
dec20 = dsp.FIRDecimator(D, h20);
dec15 = dsp.FIRDecimator(D, h15);

n = (0:Ntotal-1).';  % oscillator index

%% -------- CP lag definitions (self-contained) --------
FsWLAN_ref   = 20e6;
NfftWLAN_ref = 64;

FsLTE_ref   = 30.72e6;
NfftLTE_ref = 2048;

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

fprintf('CP lags @ FsProc=%.2f MHz: NR=%d, LTE=%d, WLAN=%d samples\n', FsProc/1e6, L_NR_p, L_LTE_p, L_WLAN_p);

numSymWin = 1;
hopFactor = 0.25;
symLenUse = 4096;
Nwin = max(256, round(numSymWin * symLenUse));
hop  = max(1,   round(hopFactor * symLenUse));

%% -------- Radar detector knobs --------
PRI_range_s = [1e-3 2e-3];
PW_range_s  = [50e-6 100e-6];
FsEnv = 2e6;

%% -------- Precompute templates at FsCommon --------
tmpl = makeTemplatesAtFsCommon(FsCommon, Ntotal, carrier);

%% -------- ROC storage (pixel samples only) --------
labelsAll = cell(C,1);
scoresAll = cell(C,1);
for c = 1:C
    labelsAll{c} = false(0,1);
    scoresAll{c} = zeros(0,1);
end

%% -------- Confusion matrix accumulator (0=Noise, 1=NR, 2=LTE, 3=WLAN, 4=RADAR) --------
classNames = {'Noise','NR','LTE','WLAN','RADAR'};
CM = zeros(5,5);          % rows=true, cols=pred
tauPix = 0;               % set >0 to reject weak stripes as Noise (e.g., 1e-4)

%% ===================== Monte Carlo loop =====================
for i = 1:Nmc
    % ---- Choose random subset of technologies ----
    present = rand(1,C) < pPresent;
    if sum(present) < minNumTech
        present(randi(C)) = true;
    end
    presentNames = classes(present);

    % ---- Random non-overlapping band allocation ----
    %[f0Map, bandsTrueHz, labelsTrue] = allocateBandsNumeric(presentNames, BW, fMin, fMax, guardHz);


    % ---- Random overlapping band allocation ----
    [f0Map, bandsTrueHz, labelsTrue] = allocateBandsAllowOverlap(presentNames, BW, fMin, fMax);

    % ---- Build xMix ----
    noise = (randn(Ntotal,1,'single') + 1j*randn(Ntotal,1,'single'))/sqrt(2);
    xMix  = noise;

    for kT = 1:numel(presentNames)
        tech = char(presentNames(kT));
        f0   = f0Map.(tech);

        x0 = single(tmpl.(tech)(:));
        x0 = x0 ./ (rms(x0)+eps('single')) * sqrt(10^(snr_dB.(tech)/10));

        osc = exp(1j*2*pi*(f0/FsCommon)*n);
        xMix = xMix + x0 .* single(osc);
    end
    xMix = single(xMix(:));

    % ---- Spectrogram ----
    [Sdb, fAxis, tAxis] = mySpectrogramDB(xMix, FsCommon, stftWin, stftHop, stftNfft);
    Nf = numel(fAxis);
    Nt = numel(tAxis);

    % ---- Ground-truth per frequency bin (bands span full time) ----
    trueOfFreq = zeros(Nf,1,'uint8');   % 0=noise, 1..4 for [NR LTE WLAN RADAR]
    for kT = 1:size(bandsTrueHz,1)
        fL = bandsTrueHz(kT,1);  fH = bandsTrueHz(kT,2);
        tech = labelsTrue(kT);
        id = classId(tech);      % must return: NR=1, LTE=2, WLAN=3, RADAR=4
        idxF = (fAxis >= fL) & (fAxis <= fH);
        trueOfFreq(idxF) = uint8(id);
    end

    % ---- Detect stripes/bands from spectrogram ----
    bandsDetHz = detectBandsFromSpectrogram(Sdb, fAxis);   % [Kdet x 2] in Hz
    Kdet = size(bandsDetHz,1);

    % ---- Stripe scoring ----
    stripeScores = zeros(Kdet, C);       % rows: stripes, cols: [NR LTE WLAN RADAR]
    stripeOfFreq = zeros(Nf,1,'uint16'); % freq-bin -> stripe index (0 if none)

    for b = 1:Kdet
        fL = bandsDetHz(b,1); fH = bandsDetHz(b,2);
        fc = 0.5*(fL+fH);
        bwEst = (fH-fL);

        idxF = (fAxis >= fL) & (fAxis <= fH);
        stripeOfFreq(idxF) = uint16(b);

        if bwEst <= 16e6
            decObj = dec15;
            FsBand = FsProc;
        else
            decObj = dec20;
            FsBand = FsProc;
        end
        reset(decObj);

        osc = exp(-1j*2*pi*(fc/FsCommon)*n);
        xBand = decObj( double(xMix) .* osc );

        rhoNR   = slidingCpCorrCoeffFast(xBand, FsBand, L_NR_p,   Nwin, hop);
        rhoLTE  = slidingCpCorrCoeffFast(xBand, FsBand, L_LTE_p,  Nwin, hop);
        rhoWLAN = slidingCpCorrCoeffFast(xBand, FsBand, L_WLAN_p, Nwin, hop);
        sRadar  = radarScoreFast(xBand, FsBand, PRI_range_s, PW_range_s, FsEnv);

        stripeScores(b,:) = [mean(rhoNR), mean(rhoLTE), mean(rhoWLAN), sRadar];
    end

    % ===== NEW (CORE FIX): predicted label is defined for the whole stripe frequency region =====
    predStripeId   = zeros(Kdet,1,'uint8');   % 1..4 (NR/LTE/WLAN/RADAR)
    bestStripeScore= zeros(Kdet,1);           % confidence per stripe

    for b = 1:Kdet
        [bestStripeScore(b), idxBest] = max(stripeScores(b,:));
        predStripeId(b) = uint8(idxBest);     % 1..4
        if bestStripeScore(b) <= tauPix
            predStripeId(b) = uint8(0);       % reject => Noise
        end
    end

    % predicted class per frequency bin (0..4)
    predOfFreq = zeros(Nf,1,'uint8');
    for jf = 1:Nf
        b = stripeOfFreq(jf);
        if b > 0
            predOfFreq(jf) = predStripeId(b);
        end
    end

    % ---- Pixel sampling ----
    % ---- Pixel sampling (STRATIFIED across 5 GT classes) ----
    Ns = min(NsampPerTrial, Nt*Nf);
    
    % classes in GT are: 0=Noise,1=NR,2=LTE,3=WLAN,4=RADAR
    ord = uint8(0:4);
    K5  = numel(ord);
    
    % target samples per class (roughly equal)
    NsPer = floor(Ns / K5);
    NsRem = Ns - NsPer*K5;
    
    it = zeros(Ns,1);   % time indices 1..Nt
    jf = zeros(Ns,1);   % freq indices 1..Nf
    
    pos = 1;
    for kk = 1:K5
        cls = ord(kk);
    
        % all frequency bins belonging to this class (bands span full time in your model)
        idxF = find(trueOfFreq == cls);
    
        % if a class is absent in this trial, skip it (reassign its budget later)
        if isempty(idxF)
            continue;
        end
    
        nTake = NsPer + (kk <= NsRem);   % distribute remainder
    
        % sample frequency bins WITH replacement (robust even if idxF is small)
        pickF = idxF(randi(numel(idxF), nTake, 1));
    
        % pick random times uniformly
        pickT = randi(Nt, nTake, 1);
    
        jf(pos:pos+nTake-1) = pickF;
        it(pos:pos+nTake-1) = pickT;
        pos = pos + nTake;
    end
    
    % If some classes were absent -> we may have underfilled; top-up from Noise (or all bins)
    if pos <= Ns
        nFill = Ns - pos + 1;
    
        idxFnoise = find(trueOfFreq == uint8(0));
        if isempty(idxFnoise)
            idxFnoise = (1:Nf).';   % fallback
        end
    
        jf(pos:Ns) = idxFnoise(randi(numel(idxFnoise), nFill, 1));
        it(pos:Ns) = randi(Nt, nFill, 1);
    end
    
    % Now sampled GT/pred labels
    gtId_u8   = trueOfFreq(jf);   % 0..4
    predId_u8 = predOfFreq(jf);   % 0..4


    % ---- Update confusion matrix ----
    %CM = CM + confusionmat(double(gtId_u8), double(predId_u8), 'Order', 0:4);
    gtAll   = double(trueOfFreq(:));   % Nf x 1
    predAll = double(predOfFreq(:));   % Nf x 1

    CM = CM + Nt * confusionmat(gtAll, predAll, 'Order', 0:4);

    if doPrintDetailTrial
        %%%%%%%Print out the statements for each trial %%%%%%%%%%%%
        % ---- Per-trial pixel counts (GT vs Pred) ----
        names5 = {'Noise','NR','LTE','WLAN','RADAR'};
        ord = uint8(0:4);
    
        gtCounts   = zeros(1,5);
        predCounts = zeros(1,5);
        for k5 = 1:5
            gtCounts(k5)   = sum(gtId_u8   == ord(k5));
            predCounts(k5) = sum(predId_u8 == ord(k5));
        end
    
        fprintf('\n--- Trial %d pixel summary (Ns=%d sampled) ---\n', i, Ns);
        for k5 = 1:5
            fprintf('%6s : GT=%8d (%.2f%%) | Pred=%8d (%.2f%%)\n', ...
                names5{k5}, gtCounts(k5), 100*gtCounts(k5)/Ns, ...
                predCounts(k5), 100*predCounts(k5)/Ns);
        end
        for techId = 1:4
            tp = sum((gtId_u8 == uint8(techId)) & (predId_u8 == uint8(techId)));
            fprintf('TP %-5s = %8d (of GT=%d)\n', names5{techId+1}, tp, gtCounts(techId+1));
        end
        fprintf('-------------------------------------------\n');
    end 
        % ---- For ROC: per-class score at sampled pixels = stripe score of that freq (else 0) ----
        % (Since bands span full time, this is consistent.)
        bIdx = stripeOfFreq(jf);  % 0..Kdet
        for c = 1:C
            s = zeros(Ns,1);
            ii = find(bIdx > 0);
            if ~isempty(ii)
                bb = double(bIdx(ii));
                s(ii) = stripeScores(bb, c);
            end
            labelsAll{c} = [labelsAll{c}; (gtId_u8 == uint8(c))];
            scoresAll{c} = [scoresAll{c}; s];
        end

    % ---- Debug plot ----
    if doPlotSpec && mod(i,plotEvery)==0
        figure(i); clf;
        imagesc(fAxis/1e6, tAxis*1e3, Sdb.'); axis xy;
        xlim([fMin fMax]/1e6); ylim([0 Ttotal_ms]);
        xlabel('Frequency (MHz)'); ylabel('Time (ms)');
        title(sprintf('Trial %d/%d | TRUE={%s} | detected stripes=%d', i, Nmc, strjoin(presentNames,','), Kdet));
        colorbar; hold on;

        % TRUE bands
        yTxt = 0.97*Ttotal_ms;
        for kT = 1:size(bandsTrueHz,1)
            fL = bandsTrueHz(kT,1); fH = bandsTrueHz(kT,2);
            xline(fL/1e6,'w--','LineWidth',1.1);
            xline(fH/1e6,'w--','LineWidth',1.1);
            text(0.5*(fL+fH)/1e6, yTxt, labelsTrue(kT), 'Color','w', 'FontWeight','bold', ...
                'HorizontalAlignment','center', 'VerticalAlignment','top');
        end

        % DETECTED bands + predicted label (stripe-wise)
        yTxt2 = 0.90*Ttotal_ms;
        for b = 1:Kdet
            fL = bandsDetHz(b,1); fH = bandsDetHz(b,2);
            xline(fL/1e6,'m-','LineWidth',1.2);
            xline(fH/1e6,'m-','LineWidth',1.2);
            pid = predStripeId(b);
            if pid==0
                lbl = "pred:Noise";
            else
                lbl = "pred:"+classes(pid);
            end
            text(0.5*(fL+fH)/1e6, yTxt2, lbl, 'Color','m', 'FontWeight','bold', ...
                'HorizontalAlignment','center', 'VerticalAlignment','top');
        end
        hold off; drawnow;
    end

    fprintf('Trial %d/%d done. TRUE={%s} | stripes=%d | sampled=%d pixels\n', ...
        i, Nmc, strjoin(presentNames,','), Kdet, Ns);
end
%% ===================== ROC (pixel-level one-vs-rest) =====================
figure; hold on; grid on;
auc = zeros(C,1);

for c = 1:C
    labels = labelsAll{c};   % logical: GT positive for class c
    scores = scoresAll{c};   % numeric: baseline score for class c

    [X,Y,~,AUC] = perfcurve(labels, scores, true);
    auc(c) = AUC;
    plot(X, Y, 'LineWidth', 1.6);
end

xlabel('False Positive Rate'); ylabel('True Positive Rate');
legend( ...
    "NR (AUC="+num2str(auc(1),3)+")", ...
    "LTE (AUC="+num2str(auc(2),3)+")", ...
    "WLAN (AUC="+num2str(auc(3),3)+")", ...
    "RADAR (AUC="+num2str(auc(4),3)+")", ...
    'Location','southeast');
title('Pixel-level ROC (sampled pixels)');



%% ===================== Confusion matrix (%) =====================
figure;
cc = confusionchart(CM, classNames, 'Normalization','row-normalized');
cc.Title = 'Confusion matrix (%)';
cc.RowSummary = 'off';
cc.ColumnSummary = 'off';
