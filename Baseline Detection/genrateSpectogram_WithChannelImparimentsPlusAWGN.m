%% ===== NR + (LTE or WLAN) + RADAR composite with channel impairments + spectrogram + GT =====
% - NR and RADAR are ALWAYS present
% - Choose either LTE OR WLAN via flag USE_ALT = "LTE" or "WLAN"
% - Adds fading + Doppler similar to MathWorks helperSpecSenseTrainingData
% - Plots: spectrogram + ground truth
clc; clear; close all;
rng(1);

%% ---------------- User flag ----------------
USE_ALT = "LTE";   % "LTE" or "WLAN"

%% ---- Class names + desired pixel IDs ----
classNames   = ["Noise" "NR" "LTE" "WLAN" "RADAR" "NRLTE" "NRWLAN" "LTEWLAN"];
pixelLabelID = uint8([0, 63, 127, 191, 6, 32, 96, 160]);   % same order as classNames
K = numel(classNames);

% Internal label indices (1..8)
Noise=1; NR=2; LTE=3; WLAN=4; RADAR=5; NRLTE=6; NRWLAN=7; LTEWLAN=8;

%% ---- Global settings ----
Fs = 100e6;                 % Hz (±50 MHz span when fmaxMHz=50)
numSubframes = 40;          % 40 ms
Tframe = numSubframes*1e-3; % seconds
SNRdB = 10;                 % desired SNR for final mixture (dB)
imSize = 256;               % 256x256
nfftSpec = 4096;            % spectrogram FFT length
guard = 1e6;                % guard from Nyquist

% Channel impairment knobs (similar to helperSpecSenseTrainingData)
useFading  = true;
DopplerMin = 0; DopplerMax = 500;       % Hz
Fc         = 5e9;                       % Hz carrier for NR/WLAN channels
DopplerHz  = DopplerMin + (DopplerMax-DopplerMin)*rand;

% Delay profile choices
NR_CDL_Profile   = "CDL-C";     % if nrCDLChannel exists
NR_DelaySpread_s = 300e-9;      % seconds
LTE_Profile      = "EPA";       % used in LTE fading fallback

N = round(Fs*Tframe);
n = (0:N-1).';

%% ==========================================================
% 1) Generate ALT technology waveform: LTE or WLAN
%% ==========================================================
if USE_ALT == "LTE"
    [altWave_rs, BWALT] = genLTE(Fs, numSubframes);
    altLabel = LTE;
else
    WLAN_BW = "CBW20";  % force 20 MHz
    [altWave_rs, BWALT] = genWLAN_fullframe_with_channel(Fs, Tframe, DopplerHz, Fc, WLAN_BW);
    altLabel = WLAN;
end
altWave_rs = fixlen(altWave_rs, N);

%% ==========================================================
% 2) NR waveform (5G Toolbox) - realistic-looking + no DM-RS warnings
%% ==========================================================
carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = 15;      % kHz
carrier.NSizeGrid = 106;             % ~19.08 MHz occupied

BWNR = carrier.NSizeGrid * 12 * carrier.SubcarrierSpacing * 1e3; % Hz
numSlots = numSubframes * carrier.SlotsPerSubframe;

% knobs for sparse activity (nice spectrogram)
slotDuty      = 0.65;
prbOccFrac    = 0.45;
symOccRange   = [6 12];
symStartRange = [0 6];
dmrsSymIdx    = 2;

nrWave = complex(zeros(0,1));
pdsch = nrPDSCHConfig;
pdsch.NumLayers   = 1;
pdsch.MappingType = "A";
pdsch.Modulation  = "16QAM";

for nslot = 0:numSlots-1
    carrier.NSlot  = mod(nslot, carrier.SlotsPerFrame);
    carrier.NFrame = floor(nslot / carrier.SlotsPerFrame);

    txGrid = nrResourceGrid(carrier, 1);

    if rand < slotDuty
        nPRB = max(6, round(prbOccFrac * carrier.NSizeGrid));
        prbs = sort(randsample(0:carrier.NSizeGrid-1, nPRB));
        pdsch.PRBSet = prbs;

        Lsym = randi(symOccRange);
        s0   = randi(symStartRange);

        if dmrsSymIdx < s0, s0 = dmrsSymIdx; end
        if dmrsSymIdx > (s0 + Lsym - 1), s0 = dmrsSymIdx - (Lsym - 1); end
        s0 = max(0, min(s0, 14 - Lsym));
        pdsch.SymbolAllocation = [s0 Lsym];

        [pdschInd, pdschInfo] = nrPDSCHIndices(carrier, pdsch);
        bits = randi([0 1], pdschInfo.G, 1);
        txGrid(pdschInd) = nrPDSCH(carrier, pdsch, bits);

        dmrsInd = nrPDSCHDMRSIndices(carrier, pdsch);
        if ~isempty(dmrsInd)
            txGrid(dmrsInd) = nrPDSCHDMRS(carrier, pdsch);
        end
    end

    w = nrOFDMModulate(carrier, txGrid, 'Nfft', 4096, 'SampleRate', Fs);
    nrWave = [nrWave; w];
end
nrWave = fixlen(nrWave, N);

%% ==========================================================
% 3) RADAR waveform (PAST style long pulse LFM) - many horizontal stripes
%% ==========================================================
randi_in = @(a,b) a + randi(b-a);

C.numBursts      = randi_in(25,30);
C.pulsesPerBurst = 1;
C.Tp             = 1e-6 * randi_in(50,100);
C.B              = 1e6  * 10;                 % 10 MHz
C.PRI            = 1e-6 * randi_in(1000,2000);
C.fs             = max(80e6, 4*C.B);

radarBase = generateType5RadarPAST(C, Tframe);
[pR,qR] = rat(Fs/C.fs, 1e-12);
radar_rs = resample(radarBase, pR, qR);
radar_rs = fixlen(radar_rs, N);
BWradar  = C.B;

%% ==========================================================
% 4) Apply channel impairments (fading + Doppler) like helperSpecSenseTrainingData
%% ==========================================================
if useFading
    nrWave_ch = applyNRChannel(nrWave, Fs, DopplerHz, Fc, NR_CDL_Profile, NR_DelaySpread_s);

    if USE_ALT == "LTE"
        altWave_ch = applyLTEChannel(altWave_rs, Fs, DopplerHz, LTE_Profile);
    else
        altWave_ch = altWave_rs; % WLAN already faded inside genWLAN_fullframe_with_channel
    end
else
    nrWave_ch  = nrWave;
    altWave_ch = altWave_rs;
end
radar_ch = radar_rs;  % keep radar clean (optional to add fading)

%% ==========================================================
% 5) Frequency placement: NR overlaps with ALT, RADAR separate
%% ==========================================================
overlapFrac = 0.30;
overlapBW   = overlapFrac * min(BWNR, BWALT);
d = (BWNR + BWALT)/2 - overlapBW;
d = max(d, 0);

maxOff_comm = Fs/2 - max(BWNR, BWALT)/2 - guard;
halfSep = min(d/2, maxOff_comm);

fOffNR  = -halfSep;
fOffALT = +halfSep;

nrShift   = nrWave_ch   .* exp(1j*2*pi*(fOffNR /Fs)*n);
altShift  = altWave_ch  .* exp(1j*2*pi*(fOffALT/Fs)*n);

% RADAR placement: fixed near +20 MHz (clamped)
fOffRADAR_des = 20e6;
maxOff_r = Fs/2 - (BWradar/2) - guard;
fOffRADAR = max(-maxOff_r, min(fOffRADAR_des, maxOff_r));
radarShift = radar_ch .* exp(1j*2*pi*(fOffRADAR/Fs)*n);

% Normalize components so mixture scaling is stable
nrShift    = nrShift   / (rms(nrShift)   + eps);
altShift   = altShift  / (rms(altShift)  + eps);
radarShift = radarShift/ (rms(radarShift)+ eps);

wNR = 1.0;
wALT = 1.0;
wR  = 0.9;

rxClean = wNR*nrShift + wALT*altShift + wR*radarShift;

%% ==========================================================
% 6) Add AWGN at SNRdB (overall mixture)
%% ==========================================================
sigP = mean(abs(rxClean).^2);
noiseP = sigP / (10^(SNRdB/10));
wgn = sqrt(noiseP/2) * (randn(N,1) + 1j*randn(N,1));
rx = rxClean + wgn;

%% ==========================================================
% 7) Plot component magnitudes
%% ==========================================================
decim = max(1, round(N/8000));
t_ms = (0:N-1).'/Fs*1e3;

% figure('Color','w');
% plot(t_ms(1:decim:end), abs(wNR*nrShift(1:decim:end))); hold on;
% plot(t_ms(1:decim:end), abs(wALT*altShift(1:decim:end)));
% plot(t_ms(1:decim:end), abs(wR*radarShift(1:decim:end)));
% grid on;
% xlabel('Time (ms)'); ylabel('|x(t)|');
% title(sprintf('Component magnitudes (decimated), Doppler=%.1f Hz, SNR=%.1f dB', DopplerHz, SNRdB));
% legend('NR', char(USE_ALT), 'RADAR', 'Location','best');

%% ==========================================================
% 8) Spectrogram image over -50..+50 MHz
%% ==========================================================
[fAxisMHz, tAxisMs, I] = makeSpecImage256(rx, Fs, nfftSpec, imSize, 50);

%% ==========================================================
% 9) Ground-truth mask only 
%% ==========================================================
gtIdx = Noise*ones(imSize, imSize, 'uint8');

fHz  = fAxisMHz(:)*1e6;
nrOn  = abs(fHz - fOffNR )  <= BWNR/2;
altOn = abs(fHz - fOffALT)  <= BWALT/2;

% RADAR mask auto-alignment (frequency)
[Srad, Frad, ~] = spectrogram(radarShift, nfftSpec, 0, nfftSpec, Fs, 'centered');
Prad = mean(abs(Srad).^2, 2);
fRadarCentHz = sum(Frad(:).*Prad(:)) / (sum(Prad(:)) + eps);
manualExtraShiftMHz = 0.0;  % if needed: 0.2..1.0
fRadarMaskCenterHz  = fRadarCentHz + manualExtraShiftMHz*1e6;
radarOn = abs(fHz - fRadarMaskCenterHz) <= BWradar/2;

for k = 1:imSize
    if radarOn(k)
        gtIdx(:,k) = RADAR;  % RADAR overrides
    else
        if nrOn(k) && altOn(k)
            gtIdx(:,k) = ternary(USE_ALT=="LTE", NRLTE, NRWLAN);
        elseif nrOn(k)
            gtIdx(:,k) = NR;
        elseif altOn(k)
            gtIdx(:,k) = altLabel; % LTE or WLAN
        end
    end
end

% Optional: map to your pixel IDs and save
gtPix = pixelLabelID(double(gtIdx));
% imwrite(gtPix, 'gtMask_pixelIDs.png');

%% ==========================================================
% 10) Plot: spectrogram + GT
%% ==========================================================
maskCmap = [ ...
    0.00 0.85 0.85;  % Noise
    0.20 0.60 1.00;  % NR
    0.55 0.40 1.00;  % LTE
    1.00 0.00 1.00;  % WLAN
    0.70 0.70 0.70;  % RADAR
    0.30 0.30 1.00;  % NRLTE
    0.90 0.60 0.00;  % NRWLAN
    0.80 0.00 0.40]; % LTEWLAN

BASE_FNT  = 15;
LABEL_FNT = 17;
TITLE_FNT = 18;
CB_FNT    = 14;

figure('Color','w','Units','inches','Position',[1 1 7.2 6.2]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
imagesc(fAxisMHz, tAxisMs, I); set(gca,'YDir','normal');
xlabel('Frequency (MHz)','FontSize',LABEL_FNT,'FontWeight','bold');
ylabel('Time (ms)','FontSize',LABEL_FNT,'FontWeight','bold');
title('Received spectrogram','FontSize',TITLE_FNT,'FontWeight','bold');
ax = gca; set(ax,'FontSize',BASE_FNT,'FontWeight','bold','LineWidth',1.2);
cb = colorbar; set(cb,'FontSize',CB_FNT,'FontWeight','bold','LineWidth',1.0);

nexttile;
imagesc(fAxisMHz, tAxisMs, gtIdx); set(gca,'YDir','normal');
xlabel('Frequency (MHz)','FontSize',LABEL_FNT,'FontWeight','bold');
ylabel('Time (ms)','FontSize',LABEL_FNT,'FontWeight','bold');
title('Ground truth','FontSize',TITLE_FNT,'FontWeight','bold');
colormap(gca, maskCmap); caxis([1 K]);
ax = gca; set(ax,'FontSize',BASE_FNT,'FontWeight','bold','LineWidth',1.2);
cb = colorbar; cb.Ticks = 1:K; cb.TickLabels = cellstr(classNames);
set(cb,'FontSize',CB_FNT,'FontWeight','bold','LineWidth',1.0);

%% ========================= Local helpers =========================
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function x = fixlen(x, N)
    x = x(:);
    if numel(x) >= N
        x = x(1:N);
    else
        x = [x; zeros(N-numel(x),1,'like',x)];
    end
end

function [altWave_rs, BWLTE] = genLTE(Fs, numSubframes)
    if exist('lteRMCDL','file')==2 && exist('lteRMCDLTool','file')==2
        rc = 'R.3';
        rmc = lteRMCDL(rc);
        rmc.TotSubframes = numSubframes;

        txDataLTE = randi([0 1], 20000, 1);
        [lteWave, ~, rmcOut] = lteRMCDLTool(rmc, txDataLTE);
        FsLTE = rmcOut.SamplingRate;
        BWLTE = rmcOut.NDLRB * 180e3;

        [p,q] = rat(Fs/FsLTE, 1e-12);
        altWave_rs = resample(lteWave, p, q);
    else
        BWLTE = 20e6;
        altWave_rs = simpleOFDMlike(Fs, 15e3, 2048, 144, 40e-3);
    end
end

function y = applyNRChannel(x, Fs, doppler, fc, cdlProfile, delaySpread_s)
    x = x(:);
    try
        if exist('nrCDLChannel','class') == 8
            chan = nrCDLChannel;
            chan.DelayProfile = char(cdlProfile);
            if isprop(chan,'DelaySpread'), chan.DelaySpread = delaySpread_s; end
            chan.SampleRate = Fs;
            chan.MaximumDopplerShift = doppler;
            chan.CarrierFrequency = fc;
            chan.Seed = randi(1e4);
            chan.TransmitAntennaArray.Size = [1 1 1 1 1];
            chan.ReceiveAntennaArray.Size  = [1 1 1 1 1];
            y = chan(x);
            return;
        end
    catch
    end

    [pd, pg] = epaTaps();
    rch = comm.RayleighChannel( ...
        'SampleRate', Fs, ...
        'PathDelays', pd, ...
        'AveragePathGains', pg, ...
        'MaximumDopplerShift', doppler, ...
        'NormalizePathGains', true);
    y = rch(x);
end

function y = applyLTEChannel(x, Fs, doppler, ~)
    x = x(:);
    try
        if exist('lteFadingChannel','file') == 2
            ch = struct();
            [pd, pg] = epaTaps();
            ch.DelayProfile = 'Custom';
            ch.PathDelays = pd;
            ch.AveragePathGaindB = pg;
            ch.NRxAnts = 1;
            ch.Seed = randi(1e4);
            ch.InitPhase = 'Random';
            ch.NTerms = 16;
            ch.InitTime = 0;
            ch.SamplingRate = Fs;
            ch.DopplerFreq = doppler;
            ch.MIMOCorrelation = 'Low';
            y = lteFadingChannel(ch, x);
            return;
        end
    catch
    end

    [pd, pg] = epaTaps();
    rch = comm.RayleighChannel( ...
        'SampleRate', Fs, ...
        'PathDelays', pd, ...
        'AveragePathGains', pg, ...
        'MaximumDopplerShift', doppler, ...
        'NormalizePathGains', true);
    y = rch(x);
end

function [pd, pg] = epaTaps()
    pd = [0 30 70 90 110 190 410]*1e-9;
    pg = [0 -1 -2 -3 -8 -17.2 -20.8];
end

function [fMHz, tMs, I] = makeSpecImage256(x, Fs, nfft, imSize, fmaxMHz)
    win = nfft;
    noverlap = 0;
    [S,F,T] = spectrogram(x, win, noverlap, nfft, Fs, 'centered');

    P = 10*log10(abs(S).^2 + eps);
    fMHz_full = F(:)/1e6;

    keep = (fMHz_full >= -fmaxMHz) & (fMHz_full <= +fmaxMHz);
    P = P(keep,:);
    fMHz_full = fMHz_full(keep);

    tMs_full = -T(:)*1e3;

    lo = prctile(P(:), 5);
    hi = prctile(P(:), 99);
    Pn = (P - lo) / max(hi-lo, eps);
    Pn = max(0, min(1, Pn));

    [Fgrid, Tgrid] = meshgrid(fMHz_full, tMs_full);
    Pfor = Pn.';

    fMHz = linspace(fMHz_full(1), fMHz_full(end), imSize);
    tMs  = linspace(tMs_full(1),  tMs_full(end),  imSize);
    [Fq, Tq] = meshgrid(fMHz, tMs);

    I = interp2(Fgrid, Tgrid, Pfor, Fq, Tq, 'linear', 0);
end

function x = generateType5RadarPAST(C, Tframe)
    fs  = C.fs;
    Tp  = C.Tp;
    B   = C.B;
    PRI = C.PRI;

    numBursts      = C.numBursts;
    pulsesPerBurst = C.pulsesPerBurst;

    Ntot = round(Tframe*fs);

    lfm = phased.LinearFMWaveform( ...
        'SampleRate', fs, ...
        'SweepBandwidth', B, ...
        'PulseWidth', Tp, ...
        'PRF', 1/PRI, ...
        'SweepDirection', 'Up', ...
        'OutputFormat', 'Pulses', ...
        'NumPulses', 1);

    pulse = lfm(); pulse = pulse(:);
    Lp = numel(pulse);

    Lpri = max(Lp+1, round(PRI*fs));
    Lburst = (pulsesPerBurst-1)*Lpri + Lp;

    oneBurst = complex(zeros(Lburst,1));
    for p = 0:pulsesPerBurst-1
        idx0 = p*Lpri + 1;
        oneBurst(idx0:idx0+Lp-1) = pulse;
    end

    x = complex(zeros(Ntot,1));
    pos = 1;
    for b = 1:numBursts
        if pos + Lburst - 1 > Ntot, break; end
        x(pos:pos+Lburst-1) = x(pos:pos+Lburst-1) + oneBurst;

        if b < numBursts
            gapSamp = round(0.001*fs*rand);
            pos = pos + Lburst + gapSamp;
        end
    end
    x = x .* exp(1j*2*pi*rand);
end

function x = simpleOFDMlike(Fs, ~, Nfft, Ncp, dur)
    N = round(dur*Fs);
    nSym = ceil(N / (Nfft+Ncp));
    x = complex(zeros(nSym*(Nfft+Ncp),1));
    for k = 1:nSym
        X = (randn(Nfft,1)+1j*randn(Nfft,1))/sqrt(2);
        s = ifft(X);
        s = [s(end-Ncp+1:end); s];
        x((k-1)*(Nfft+Ncp)+1:k*(Nfft+Ncp)) = s;
    end
    x = x(1:N);
end

function [wlan_rs, BWWLAN] = genWLAN_fullframe_with_channel(FsCommon, Tframe, dopplerHz, fcHz, chBW)
    cfg = wlanNonHTConfig('ChannelBandwidth', char(chBW), 'MCS', 4, 'PSDULength', 1000);
    bits = randi([0 1], cfg.PSDULength*8, 1);
    tx   = wlanWaveformGenerator(bits, cfg);
    FsW  = wlanSampleRate(cfg);

    Nnative = round(Tframe * FsW);
    txFull  = repToLen(tx, Nnative);

    rxFull = applyWLANChannel_native(txFull, FsW, dopplerHz, fcHz, chBW);

    [p,q] = rat(FsCommon/FsW, 1e-12);
    wlan_rs = resample(rxFull, p, q);

    if chBW == "CBW20"
        BWWLAN = 20e6;
    else
        BWWLAN = 40e6;
    end
end

function y = applyWLANChannel_native(x, FsW, dopplerHz, fcHz, chBW)
    x = x(:);

    try
        if exist('wlanTGaxChannel','file') == 2
            ch = wlanTGaxChannel;
            if isprop(ch,'SampleRate'),          ch.SampleRate = FsW; end
            if isprop(ch,'CarrierFrequency'),    ch.CarrierFrequency = fcHz; end
            if isprop(ch,'MaximumDopplerShift'), ch.MaximumDopplerShift = dopplerHz; end
            if isprop(ch,'DelayProfile'),        ch.DelayProfile = 'Model-B'; end
            if isprop(ch,'ChannelBandwidth'),    ch.ChannelBandwidth = char(chBW); end
            if isprop(ch,'Seed'),                ch.Seed = randi(1e4); end
            y = ch(x);
            return;
        end
    catch
    end

    [pd, pg] = epaTaps();
    rch = comm.RayleighChannel( ...
        'SampleRate', FsW, ...
        'PathDelays', pd, ...
        'AveragePathGains', pg, ...
        'MaximumDopplerShift', dopplerHz, ...
        'NormalizePathGains', true);
    y = rch(x);
end

function y = repToLen(x, N)
    x = x(:);
    if isempty(x)
        y = zeros(N,1,'like',x);
        return;
    end
    R = ceil(N / numel(x));
    y = repmat(x, R, 1);
    y = y(1:N);
end