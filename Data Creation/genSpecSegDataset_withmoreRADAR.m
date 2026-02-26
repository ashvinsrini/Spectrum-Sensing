%% ===================== DATASET: SAVE COLORED SPECTROGRAM + LABEL MASK =====================
clc; close all; clear;
rng(1);

%% -------- Dataset knobs --------
Nmc = 1;                       % number of spectrogram/mask pairs to generate
radarPulseWiden_us = 200;   % <<< widen each radar pulse by +/-200 us (tune: 80..300)

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
labelVals = uint8([0 63 127 191 6]);   % index = trueOfFreq+1

% Enhanced mask visualization colors (ORDER MUST MATCH [Noise NR LTE WLAN RADAR])
maskColors_u8 = uint8(round(255 * [ ...
    0.20 0.90 0.90;   % Noise
    0.10 0.60 1.00;   % NR
    0.70 0.40 0.90;   % LTE
    1.00 0.00 1.00;   % WLAN
    1.00 0.65 0.00;   % RADAR
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
numRadarRegions = 1;      % <<< REQUIRED: 3 vertical RADAR regions every spectrogram
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

%% ===================== Monte Carlo loop =====================
for i = 1:Nmc

    % ==========================================================
    % 1) Choose OTHER technologies (NR/LTE/WLAN) uniformly.
    %    RADAR will be enforced as 3 separate non-overlapping bands.
    % ==========================================================
    otherClasses = ["NR","LTE","WLAN"];
    Kother = randi([0 numel(otherClasses)], 1, 1);     % allow 0..3 other techs
    idxO   = randperm(numel(otherClasses), Kother);
    presentOther = otherClasses(idxO);
    
    % ==========================================================
    % 2) Allocate NON-OVERLAPPING bands for:
    %    presentOther + (RADAR repeated 3 times)
    %    Returns:
    %      labelsTrue  : Kbands x 1 string array (RADAR repeats)
    %      bandsTrueHz : Kbands x 2 [fL fH]
    %      f0ListHz    : Kbands x 1 center freqs
    % ==========================================================
    % [labelsTrue, bandsTrueHz, f0ListHz] = allocateBandsNonOverlap_multiRadar_local( ...
    %     presentOther, BW, fMin, fMax, guardHz, numRadarRegions);
    [labelsTrue, bandsTrueHz, f0ListHz] = allocateBandsRadarExclusive_local( ...
    presentOther, BW, fMin, fMax, guardHz, numRadarRegions, allowOverlapNonRadar);

    presentNamesUnique = unique([presentOther "RADAR"], 'stable');  % for printing only
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

    % 5) Ground-truth label per freq bin (0..4 using your classId)
    trueOfFreq = zeros(Nf,1,'uint8');   % 0=noise, 1..4 for [NR LTE WLAN RADAR]
    for kT = 1:size(bandsTrueHz,1)
        fL = bandsTrueHz(kT,1);  fH = bandsTrueHz(kT,2);
        tech = labelsTrue(kT);
        id = classId(tech);      % NR=1, LTE=2, WLAN=3, RADAR=4
        idxF = (fAxis >= fL) & (fAxis <= fH);
        trueOfFreq(idxF) = uint8(id);
    end

    % 6) Build MASK (Nt x Nf) with REQUIRED values {0,63,127,191,6}
    % Start with the "frequency-only" mask (works for NR/LTE/WLAN)
    maskOfFreq = labelVals(double(trueOfFreq) + 1);   % Nf x 1
    maskNtNf   = repmat(maskOfFreq(:).', Nt, 1);      % Nt x Nf
    
    % --- FIX: make RADAR time-localized using TF energy inside radar band ---
    % This will ONLY label radar pixels where spectrogram energy is high;
    % remaining pixels inside the radar band become Noise.
    maskNtNf = refineRadarMaskTF_local(maskNtNf, Sdb, trueOfFreq, labelVals);

    % 7) Colored spectrogram PNG like your debug plot (uses Sdb.')
    specRGB = SdbToRGB_local(Sdb, cmap, useFixedCLim, clim_dB); % Nt x Nf x 3 (uint8)

    % Match axis xy look
    specRGB  = flipud(specRGB);
    maskNtNf = flipud(maskNtNf);

    % Resize
    if ~isempty(outSize)
        specRGB  = imresize(specRGB,  outSize, 'bilinear');
        maskNtNf = imresize(maskNtNf, outSize, 'nearest');
    end

    % Save PNG spectrogram + PNG mask
    imwrite(specRGB,          fullfile(specDir, sprintf('spec_%04d.png', i)));
    imwrite(uint8(maskNtNf),  fullfile(maskDir, sprintf('mask_%04d.png', i)));

    % Save HDF mask (same uint8 values)
    hdfPath = fullfile(maskHdfDir, sprintf('mask_%04d.hdf', i));
    if exist(hdfPath,'file'), delete(hdfPath); end
    imwrite(uint8(maskNtNf), hdfPath, 'hdf');

    % Save enhanced (RGB) mask for visualization (maps by VALUES, not 1..5)
    maskRGB = maskToRGB_values_local(uint8(maskNtNf), labelVals, maskColors_u8);
    imwrite(maskRGB, fullfile(maskEnhDir, sprintf('maskenh_%04d.png', i)));

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
        sprintf('spec_%04d.png',i), sprintf('mask_%04d.png',i), sprintf('mask_%04d.hdf',i));

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