%% ===================== DATASET: SAVE COLORED SPECTROGRAM + LABEL MASK =====================
clc; close all; clear;
rng(1);

%% -------- Dataset knobs --------
Nmc = 30;                       % number of spectrogram/mask pairs to generate
radarPulseWiden_us = 200;   % <<< widen each radar pulse by +/-200 us (tune: 80..300)
pOverlap = 0.75;          % probability of overlap (you already asked this)
minOvFrac = 0.15;         % minimum overlap fraction of min(BW1,BW2)
maxOvFrac = 0.60;         % maximum overlap fraction of min(BW1,BW2)  <-- key for "not complete"

% If true: ONLY generate & save dataset (skip baseline detection/ROC/CM)
doDatasetOnly = false;
allowOverlapNonRadar = true;   % NR/LTE/WLAN may overlap each other (RADAR stays exclusive)

% Output folders
outRoot = fullfile(pwd,'dataset_out');
specDir = fullfile(outRoot,'spectrogram');   ensureDir_local(specDir);
maskDir = fullfile(outRoot,'mask');          ensureDir_local(maskDir);
maskEnhDir = fullfile(outRoot,'maskenhanced'); ensureDir_local(maskEnhDir);
maskHdfDir = fullfile(outRoot,'maskhdf');    ensureDir_local(maskHdfDir);

% Save options
outSize = [256 256];             % [] to keep native size (Nt x Nf)
useFixedCLim = true;             % fix color scaling for consistent appearance
clim_dB = [-40 60];              % adjust to taste (like your plot's range)
cmap = parula(256);              % colormap for saving RGB PNG

% ---- Desired output label values (uint8) ----
% Internal trueOfFreq convention: 0=Noise, 1=NR, 2=LTE, 3=WLAN, 4=RADAR
% Required exported values: Noise=0, NR=63, LTE=127, WLAN=191, RADAR=6
% ---- Mask label values (uint8) ----
% Base: Noise=0, NR=63, LTE=127, WLAN=191
% Overlap classes (choose any unique uint8 values):
val.NOISE   = uint8(0);
val.NR      = uint8(63);
val.LTE     = uint8(127);
val.WLAN    = uint8(191);
val.NRLTE   = uint8(32);
val.NRWLAN  = uint8(96);
val.LTEWLAN = uint8(160);

% For visualization mapping (ORDER matters here)
labelValsViz = uint8([ ...
    val.NOISE, val.NR, val.LTE, val.WLAN, val.NRLTE, val.NRWLAN, val.LTEWLAN ]);

% Enhanced mask visualization colors (7 rows)
maskColors_u8 = uint8(round(255 * [ ...
    0.20 0.90 0.90;   % Noise
    0.10 0.60 1.00;   % NR
    0.70 0.40 0.90;   % LTE
    1.00 0.00 1.00;   % WLAN
    0.10 0.90 0.10;   % NRLTE   (pick any colors you like)
    1.00 0.20 0.20;   % NRWLAN
    1.00 1.00 0.20;   % LTEWLAN
]));

noiseGain_dB = 8;     % try 6..15 dB
noiseGain    = 10^(noiseGain_dB/20);


%% -------- Monte Carlo knobs (kept for compatibility) --------
classes = ["NR","LTE","WLAN","RADAR"];
C = numel(classes);

pPresent   = 0.60;              
minNumTech = 1;                 
NsampPerTrial = 200000;         

%% -------- Wideband capture --------
FsCommon  = 300e6;
Ttotal_ms = 40;
Ntotal    = round(Ttotal_ms*1e-3*FsCommon);

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

% ---- Force MULTI-RADAR regions per image ----
numRadarRegions = 4;      % <<< REQUIRED: 3 vertical RADAR regions every spectrogram
radarRandShift  = true;   % optional: make each region's pulse pattern look different


%% -------- Spectrogram params --------
doPrintDetailTrial = false;
doPlotSpec = false;
useNonOverlap = false;
plotEvery  = 1;

stftWin  = 4096;
stftHop  = 4096;
stftNfft = 4096;

%% -------- FAST per-stripe processing: DDC + decimate --------
% (kept as-is; unused if you don't run baseline)
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
% (kept as-is; unused if you don't run baseline)
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

% ===== Two-nonRADAR-tech schedule (3 combos) + force 50% overlap =====
pairList = ["NR","LTE";
            "NR","WLAN";
            "LTE","WLAN"];     % 3 combos

Npairs = size(pairList,1);

% balance combos as evenly as possible over Nmc
base = floor(Nmc/Npairs);
rem  = Nmc - base*Npairs;

counts = base*ones(1,Npairs);
if rem > 0
    extra = randperm(Npairs, rem);
    counts(extra) = counts(extra) + 1;
end

pairPerTrial = strings(Nmc,2);
idx = 1;
for p = 1:Npairs
    if counts(p) > 0
        pairPerTrial(idx:idx+counts(p)-1,:) = repmat(pairList(p,:), counts(p), 1);
        idx = idx + counts(p);
    end
end

% exactly 50% trials overlap (for odd Nmc => closest possible)
Nover = round(pOverlap*Nmc);
overlapPerTrial = false(Nmc,1);
overlapPerTrial(1:Nover) = true;

% shuffle both together
perm = randperm(Nmc);
pairPerTrial = pairPerTrial(perm,:);
overlapPerTrial = overlapPerTrial(perm);



%% ===================== Monte Carlo loop =====================
for i = 1:Nmc

    % ==========================================================
    % ==========================================================
    % 1) Pick TWO non-RADAR technologies for this spectrogram
    % ==========================================================
    techA = string(pairPerTrial(i,1));
    techB = string(pairPerTrial(i,2));
    doOverlap = overlapPerTrial(i);   % true for 50% of trials
    
    labelsTrue  = [techA; techB];
    bandsTrueHz = zeros(2,2);
    f0ListHz    = zeros(2,1);
    
    % ==========================================================
    % 2) Place techA band anywhere
    % ==========================================================
    BW1    = BW.(techA);
    f0min1 = fMin + BW1/2 + guardHz;
    f0max1 = fMax - BW1/2 - guardHz;
    if f0min1 >= f0max1
        error('No room to place %s band.', techA);
    end
    
    f01 = f0min1 + (f0max1 - f0min1)*rand();
    fL1 = f01 - BW1/2;
    fH1 = f01 + BW1/2;
    
    bandsTrueHz(1,:) = [fL1 fH1];
    f0ListHz(1)      = f01;
    
    % ==========================================================
    % 3) Place techB band (forced overlap OR forced non-overlap)
    % ==========================================================
    BW2    = BW.(techB);
    f0min2 = fMin + BW2/2 + guardHz;
    f0max2 = fMax - BW2/2 - guardHz;
    if f0min2 >= f0max2
        error('No room to place %s band.', techB);
    end
    
    placed = false;
    
    if doOverlap
        ovMin = max(f0min2, fL1 - BW2/2);
        ovMax = min(f0max2, fH1 + BW2/2);
    
        placed = false;
        minBW  = min(BW1, BW2);
    
        for attempt = 1:800
            if ~(ovMin < ovMax), break; end
    
            f02 = ovMin + (ovMax - ovMin)*rand();
            fL2 = f02 - BW2/2;
            fH2 = f02 + BW2/2;
    
            overlapBW = min(fH1,fH2) - max(fL1,fL2);   % Hz
            if overlapBW <= 0
                continue;
            end
    
            % reject "complete overlap / containment"
            fullContain12 = (fL2 >= fL1) && (fH2 <= fH1);   % band2 inside band1
            fullContain21 = (fL1 >= fL2) && (fH1 <= fH2);   % band1 inside band2
            if fullContain12 || fullContain21
                continue;
            end
    
            % enforce "partial overlap amount"
            if overlapBW < minOvFrac*minBW, continue; end
            if overlapBW > maxOvFrac*minBW, continue; end
    
            placed = true;
            break;
        end
    
        if ~placed
            % fallback: treat as non-overlap trial if couldn't satisfy constraints
            doOverlap = false;
        end
    end

    if ~doOverlap
        for attempt = 1:800
            f02 = f0min2 + (f0max2 - f0min2)*rand();
            fL2 = f02 - BW2/2;
            fH2 = f02 + BW2/2;
    
            if (fH2 < fL1 - guardHz) || (fL2 > fH1 + guardHz)
                placed = true;
                break;
            end
        end
    end
    
    if ~placed
        error('Could not place %s relative to %s (doOverlap=%d). Reduce guardHz.', techB, techA, doOverlap);
    end
    
    bandsTrueHz(2,:) = [fL2 fH2];
    f0ListHz(2)      = f02;
    
    presentNamesUnique = unique(labelsTrue, 'stable');
    presentNames = presentNamesUnique;


    % ==========================================================
    % 3) Build xMix by summing each allocated band (RADAR appears 3 times)
    % ==========================================================
    %noise = noiseGain * (randn(Ntotal,1,'single') + 1j*randn(Ntotal,1,'single'))/sqrt(2);
    noise = 0.001*(randn(Ntotal,1,'single') + 1j*randn(Ntotal,1,'single'))/sqrt(2);
    xMix  = noise;

    
    for kBand = 1:numel(labelsTrue)
        tech = char(labelsTrue(kBand));
        f0   = f0ListHz(kBand);
    
        x0 = single(tmpl.(tech)(:));
    
        % Optional: make the 3 radar regions have different pulse timing
        % (keeps "multiple horizontal stripes" but avoids identical copies)
        if strcmpi(tech,'RADAR') && radarRandShift
            sh = randi([0, Ntotal-1], 1, 1);
            x0 = circshift(x0, sh);
            x0 = x0 .* exp(1j*2*pi*rand());   % random phase
        end
        % --- NEW: widen RADAR horizontal stripes in time domain ---
        if strcmpi(tech,'RADAR')
            x0 = widenRadarPulses_local(x0, FsCommon, radarPulseWiden_us);
        end
        % Scale to desired SNR
        x0 = x0 ./ (rms(x0)+eps('single')) * sqrt(10^(snr_dB.(tech)/10));
    
        osc  = exp(1j*2*pi*(f0/FsCommon)*n);
        xMix = xMix + x0 .* single(osc);
    end
    xMix = single(xMix(:));


    % 4) Spectrogram (same as your pipeline)
    [Sdb, fAxis, tAxis] = mySpectrogramDB(xMix, FsCommon, stftWin, stftHop, stftNfft);
    Nf = numel(fAxis);
    Nt = numel(tAxis);

    % ==========================================================
    % 5+6) Build masks with rule:
    %      - partial overlap bins are assigned to LEFT-most tech (lower fL)
    %      - produces 2 solid regions in PNG, and multi-hot is "resolved" too
    % ==========================================================
    % ==========================================================
    % 5+6) Build masks with overlap encoded as NEW overlap-classes
    %      - non-overlap pixels keep base tech labels
    %      - overlap pixels become one of: NRLTE / NRWLAN / LTEWLAN
    % ==========================================================
    techNames = ["NR","LTE","WLAN"];      % 3 channels
    Cmask = numel(techNames);
    
    % bands
    fL1 = bandsTrueHz(1,1);  fH1 = bandsTrueHz(1,2);
    fL2 = bandsTrueHz(2,1);  fH2 = bandsTrueHz(2,2);
    
    tech1 = string(labelsTrue(1));
    tech2 = string(labelsTrue(2));
    
    idx1  = (fAxis >= fL1) & (fAxis <= fH1);
    idx2  = (fAxis >= fL2) & (fAxis <= fH2);
    idxOv = idx1 & idx2;
    
    idx1o = idx1 & ~idxOv;   % only tech1
    idx2o = idx2 & ~idxOv;   % only tech2
    
    % ---- Multi-hot mask (TRUE overlap => both channels true in overlap) ----
    trueMaskF = false(Nf, Cmask);  % Nf x 3
    ch1 = find(techNames == tech1, 1);
    ch2 = find(techNames == tech2, 1);
    if isempty(ch1) || isempty(ch2)
        error('Unexpected tech(s): %s and %s', tech1, tech2);
    end
    trueMaskF(idx1, ch1) = true;   % includes overlap
    trueMaskF(idx2, ch2) = true;   % includes overlap
    
    maskMH = repmat(reshape(trueMaskF, 1, Nf, Cmask), Nt, 1, 1);
    
    % ---- Single-channel PNG mask with overlap as combo-class ----
    maskOfFreq_png = zeros(1, Nf, 'uint8');
    
    % base tech mapping
    maskOfFreq_png(idx1o) = val.(char(tech1));   % NR/LTE/WLAN -> 63/127/191
    maskOfFreq_png(idx2o) = val.(char(tech2));
    
    % overlap mapping
    if ( (tech1=="NR"   && tech2=="LTE")  || (tech1=="LTE"  && tech2=="NR") )
        ovVal = val.NRLTE;
    elseif ( (tech1=="NR"  && tech2=="WLAN") || (tech1=="WLAN" && tech2=="NR") )
        ovVal = val.NRWLAN;
    elseif ( (tech1=="LTE" && tech2=="WLAN") || (tech1=="WLAN" && tech2=="LTE") )
        ovVal = val.LTEWLAN;
    else
        error('Unexpected overlap pair: %s + %s', tech1, tech2);
    end
    
    maskOfFreq_png(idxOv) = ovVal;
    
    maskNtNf = repmat(maskOfFreq_png, Nt, 1); % Nt x Nf

    % 7) Colored spectrogram PNG like your debug plot (uses Sdb.')
    specRGB = SdbToRGB_local(Sdb, cmap, useFixedCLim, clim_dB); % Nt x Nf x 3 (uint8)

    % Match axis xy look
    specRGB  = flipud(specRGB);
    maskNtNf = flipud(maskNtNf);
    maskMH = flipud(maskMH);


    % Resize
    if ~isempty(outSize)
        specRGB  = imresize(specRGB,  outSize, 'bilinear');
        maskNtNf = imresize(maskNtNf, outSize, 'nearest');
        % resize multi-hot mask per channel
        maskMH = imresize(uint8(maskMH), outSize, 'nearest');
        maskMH = logical(maskMH);

    end

    % Save PNG spectrogram + PNG mask
    imwrite(specRGB,          fullfile(specDir, sprintf('spec_twotechNonRADAR%04d.png', i)));
    imwrite(uint8(maskNtNf),  fullfile(maskDir, sprintf('mask_twotechNonRADAR%04d.png', i)));

    % Save HDF mask (same uint8 values)
    hdfPath = fullfile(maskHdfDir, sprintf('mask_twotechNonRADAR%04d.hdf', i));
    if exist(hdfPath,'file'), delete(hdfPath); end
    imwrite(uint8(maskNtNf), hdfPath, 'hdf');

    % Save enhanced (RGB) mask for visualization (maps by VALUES, not 1..5)
    maskRGB = maskToRGB_values_local(uint8(maskNtNf), labelValsViz, maskColors_u8);
    imwrite(maskRGB, fullfile(maskEnhDir, sprintf('maskenh_twotechNonRADAR%04d.png', i)));

    % Optional debug plot
    if doPlotSpec && mod(i,plotEvery)==0
        figure(i); clf;
        imagesc(fAxis/1e6, tAxis*1e3, Sdb.'); axis xy;
        xlim([fMin fMax]/1e6); ylim([0 Ttotal_ms]);
        xlabel('Frequency (MHz)'); ylabel('Time (ms)');
        title(sprintf('Trial %d/%d | TRUE={%s}', i, Nmc, strjoin(string(presentNamesUnique),',')));
        colorbar; hold on;
        yTxt = 0.97*Ttotal_ms;
        for kT = 1:size(bandsTrueHz,1)
            fL = bandsTrueHz(kT,1); fH = bandsTrueHz(kT,2);
            xline(fL/1e6,'w--','LineWidth',1.1);
            xline(fH/1e6,'w--','LineWidth',1.1);
            text(0.5*(fL+fH)/1e6, yTxt, labelsTrue(kT), 'Color','w', 'FontWeight','bold', ...
                'HorizontalAlignment','center', 'VerticalAlignment','top');
        end
        hold off; drawnow;
    end

    fprintf('Saved %d/%d | TRUE={%s} | spec=%s | mask=%s | hdf=%s\n', ...
        i, Nmc, strjoin(presentNames,','), ...
        sprintf('spec_twotechNonRADAR%04d.png',i), sprintf('mask_twotechNonRADAR%04d.png',i), sprintf('mask_twotechNonRADAR%04d.hdf',i));

    if doDatasetOnly
        continue;
    end

    % If you keep doDatasetOnly=false and want baseline ROC/CM,
    % paste the original baseline code below unchanged.
end

fprintf('\nDONE.\nSpectrograms: %s\nMasks:        %s\nMaskEnhanced: %s\nMaskHDF:      %s\n', ...
    specDir, maskDir, maskEnhDir, maskHdfDir);



%%%%%%%%%%%%% saving spectrogram and maskhdf files in one folder 
root      = outRoot;
specDir   = fullfile(root,'spectrogram');
maskHdfDir= fullfile(root,'maskhdf');
outDir    = fullfile(root,'inputfiles');

if ~exist(outDir,'dir'); mkdir(outDir); end

% --- Copy spectrograms ---
specFiles = dir(fullfile(specDir,'spec_*.png'));
if isempty(specFiles)
    error('No spec_*.png found in %s', specDir);
end

for k = 1:numel(specFiles)
    src = fullfile(specDir, specFiles(k).name);
    dst = fullfile(outDir,  specFiles(k).name);
    copyfile(src, dst);
end

% --- Copy HDF masks ---
maskFiles = dir(fullfile(maskHdfDir,'mask_*.hdf'));
if isempty(maskFiles)
    error('No mask_*.hdf found in %s', maskHdfDir);
end

for k = 1:numel(maskFiles)
    src = fullfile(maskHdfDir, maskFiles(k).name);
    dst = fullfile(outDir,     maskFiles(k).name);
    copyfile(src, dst);
end

fprintf('Copied %d specs and %d masks into: %s\n', numel(specFiles), numel(maskFiles), outDir);