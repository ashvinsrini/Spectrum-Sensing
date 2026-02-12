%% ==============================================================
%  Composite wideband xMix: NR / LTE / WLAN in different sub-bands
%  (ALL technologies present for the FULL time 0..40 ms)
%  Goal: Spectrogram with Frequency on x-axis [-150..150] MHz, Time on y-axis [0..40] ms
%  Classes (mask): 0=Noise, 1=NR, 2=LTE, 3=WLAN
%% ==============================================================

clear; clc; close all;
rng(7);

%% -------------------- User knobs --------------------
FsCommon   = 300e6;         % => spectrogram frequency axis approx [-150,150] MHz
Ttotal_ms  = 40;            % total duration (ms)

% Desired SNR of each technology relative to unit-variance noise (over full 40 ms)
snrNR_dB   = 0;
snrLTE_dB  = 0;
snrWLAN_dB = 0;

% Sub-band centers (Hz) within [-150..150] MHz
f0_NR   = -40e6;
f0_LTE  = 10e6;
f0_WLAN = +90e6;

% Approx occupied bandwidths used for label mask only (Hz)
BW_NR   = 20e6;
BW_LTE  = 20e6;
BW_WLAN = 20e6;

% Taper only at start/end of the FULL 40 ms frame (reduces splatter at boundaries)
fade_ms = 0.2;

% Spectrogram parameters
stftWin  = 4096;
stftHop  = 4096;
stftNfft = 4096;

%% -------------------- Derived lengths --------------------
Ntotal = round(Ttotal_ms*1e-3 * FsCommon);
Nfade  = round(fade_ms*1e-3 * FsCommon);

%% =========================================================
%  1) Generate baseband waveforms at native sampling rates
% =========================================================

% ----- WLAN (802.11a-like Non-HT, 20 MHz) -----
cfgWLAN = wlanNonHTConfig('ChannelBandwidth','CBW20','MCS',4,'PSDULength',1000);
bitsWiFi = randi([0 1], cfgWLAN.PSDULength*8, 1);
txWLAN   = wlanWaveformGenerator(bitsWiFi, cfgWLAN);
FsWLAN   = wlanSampleRate(cfgWLAN);

% ----- NR OFDM (~20 MHz by using 106 PRBs @ 15 kHz) -----
carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = 30;      % kHz
carrier.NSizeGrid         = 106;     % ~19.08 MHz (106*12*15k)
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

% ----- LTE (use LTE toolbox if present, else fallback LTE-like OFDM) -----
if exist('lteRMCDL','file') == 2 && exist('lteRMCDLTool','file') == 2
    rmccfg = lteRMCDL('R.7');     % 20 MHz LTE DL RMC
    rmccfg.NCellID = 0;
    rmccfg.NSubframe = 0;
    [txLTE, ~, infoLTE] = lteRMCDLTool(rmccfg, randi([0 1], rmccfg.PDSCH.TrBlkSizes(1), 1));
    FsLTE = infoLTE.SamplingRate;
else
    FsLTE = 30.72e6;             % LTE 20 MHz nominal sampling rate
    txLTE = simpleLTElikeOFDM(FsLTE);
end

%% =========================================================
%  2) Resample all to FsCommon and tile to FULL duration (40 ms)
% =========================================================
xWLAN = resampleAny(txWLAN, FsWLAN, FsCommon);
xNR   = resampleAny(txNR,   FsNR,   FsCommon);
xLTE  = resampleAny(txLTE,  FsLTE,  FsCommon);

% Make each exactly Ntotal by repeating + truncating (ALL active for full 40 ms)
xWLAN = repToLen(xWLAN, Ntotal);
xNR   = repToLen(xNR,   Ntotal);
xLTE  = repToLen(xLTE,  Ntotal);

xWLAN = single(xWLAN);
xNR   = single(xNR);
xLTE  = single(xLTE);

%% =========================================================
%  3) Set relative SNRs (w.r.t. unit-variance wideband noise)
% =========================================================
% Wideband complex noise across full duration
noise = (randn(Ntotal,1,'single') + 1i*randn(Ntotal,1,'single'))/sqrt(2);

% Normalize each signal to unit RMS then scale to target SNR w.r.t. noise power (=1)
xNR   = xNR   ./ (rms(xNR)   + eps('single')) * sqrt(10^(snrNR_dB/10));
xLTE  = xLTE  ./ (rms(xLTE)  + eps('single')) * sqrt(10^(snrLTE_dB/10));
xWLAN = xWLAN ./ (rms(xWLAN) + eps('single')) * sqrt(10^(snrWLAN_dB/10));

% Optional taper at the frame edges only
xNR   = applyTaper(xNR,   Nfade);
xLTE  = applyTaper(xLTE,  Nfade);
xWLAN = applyTaper(xWLAN, Nfade);

%% =========================================================
%  4) Frequency-shift each technology into its own band (FULL duration)
% =========================================================
n = single((0:Ntotal-1).');

xNR_shift   = xNR   .* exp(1i*2*pi*(f0_NR/FsCommon)*n);
xLTE_shift  = xLTE  .* exp(1i*2*pi*(f0_LTE/FsCommon)*n);
xWLAN_shift = xWLAN .* exp(1i*2*pi*(f0_WLAN/FsCommon)*n);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RADAR generation using PAST Tool box

%% =========================================================
%  RADAR: Long-pulse burst LFM waveform at FsCommon, duration Ntotal
%  Produces xRADAR_shift (same length as xMix), centered at f0_RADAR
%  Requires Phased Array System Toolbox for phased.LinearFMWaveform (optional).
%  If not available, it falls back to analytic complex LFM generation.
%% =========================================================

% ---- User knobs (choose within your Table-6 ranges) ----
snrRADAR_dB    = 10;          % target SNR during pulse-active samples (w.r.t unit noise var=1)
f0_RADAR       = +90e6;       % radar band center (Hz) keep within [-150,150] MHz and avoid overlaps
pulseWidth_us  = 80;          % 50..100 us
chirpBW_MHz    = 15;          % 5..20 MHz
PRI_us         = 1500;        % 1000..2000 us
pulsesPerBurst = 2;           % 1..3
numBursts      = 12;          % 8..20

% ---- Convert to seconds / samples ----
PW   = pulseWidth_us * 1e-6;
BW   = chirpBW_MHz   * 1e6;
PRI  = PRI_us        * 1e-6;

PW_samp  = max(1, round(PW  * FsCommon));
PRI_samp = max(PW_samp+1, round(PRI * FsCommon));

% Total duration
Ttotal = Ntotal / FsCommon;

% Burst duration (approx)
burstDur = (pulsesPerBurst-1)*PRI + PW;

% Evenly place bursts across the 40 ms so all fit
% (leave small guard so last pulse doesn't exceed frame)
tStartMax = max(0, Ttotal - burstDur - 1/FsCommon);
burstStarts = linspace(0, tStartMax, numBursts);

% ---- Preallocate baseband radar waveform (centered at 0 Hz) ----
xRADAR_bb = zeros(Ntotal,1,'single');
active = false(Ntotal,1);  % pulse-active sample mask (for SNR scaling)

% ---- Build one LFM pulse generator (toolbox if possible, else analytic) ----
usePhased = (exist('phased.LinearFMWaveform','class') == 8);

if usePhased
    % Note: this object typically outputs 1 PRI worth of samples per call.
    wfm = phased.LinearFMWaveform( ...
        'SampleRate', FsCommon, ...
        'PulseWidth', PW, ...
        'PRF', 1/PRI, ...
        'SweepBandwidth', BW, ...
        'SweepDirection', 'Up');
end

% ---- Insert bursts/pulses into xRADAR_bb ----
for b = 1:numBursts
    t0 = burstStarts(b);

    for p = 1:pulsesPerBurst
        idx0 = round((t0 + (p-1)*PRI) * FsCommon) + 1;
        idx1 = idx0 + PW_samp - 1;
        if idx0 < 1 || idx1 > Ntotal
            continue; % skip anything that would exceed frame
        end

        % Generate one complex LFM pulse (length PW_samp)
        if usePhased
            % Get one PRI from system object, then take only the active pulse part
            xPRI = wfm(); xPRI = xPRI(:);
            if numel(xPRI) < PW_samp
                % fallback if object returned only PW samples
                pulse = xPRI;
            else
                pulse = xPRI(1:PW_samp);
            end
        else
            % Analytic complex baseband LFM from -BW/2 to +BW/2 over PW
            t = (0:PW_samp-1).' / FsCommon;
            k = BW / PW;                  % sweep rate (Hz/s)
            f0 = -BW/2;                   % start frequency (Hz)
            pulse = exp(1j*2*pi*( f0*t + 0.5*k*t.^2 ));
        end

        % Optional taper inside the pulse to reduce splatter (nice looking spectrogram)
        w = single(hann(numel(pulse)));
        pulse = single(pulse) .* w;

        % Add pulse into frame
        xRADAR_bb(idx0:idx1) = xRADAR_bb(idx0:idx1) + pulse;
        active(idx0:idx1) = true;
    end
end

% ---- Scale radar to desired SNR DURING active pulse samples (noise var=1) ----
if any(active)
    rmsActive = rms(double(xRADAR_bb(active)));
    xRADAR_bb = xRADAR_bb ./ (rmsActive + eps('single')) * sqrt(10^(snrRADAR_dB/10));
end

% ---- Shift radar into its sub-band ----
n = single((0:Ntotal-1).');
xRADAR_shift = xRADAR_bb .* exp(1i*2*pi*(f0_RADAR/FsCommon)*n);









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% =========================================================
%  5) Composite xMix: Noise + NRband + LTEband + WLANband (0..40 ms)
% =========================================================
%xMix = noise + xNR_shift + xLTE_shift + xWLAN_shift;
xMix = noise + xNR_shift + xLTE_shift+ xRADAR_shift;
%% =========================================================
%  6) Spectrogram (freq on x, time on y) and label mask
% =========================================================
[Sdb, f, t] = mySpectrogramDB(xMix, FsCommon, stftWin, stftHop, stftNfft);

figure;
imagesc(f/1e6, t*1e3, Sdb.'); axis xy;
xlim([-150 150]); ylim([0 Ttotal_ms]);
xlabel('Frequency (MHz)'); ylabel('Time (ms)');
title('Spectrogram (dB): Freq on x-axis, Time on y-axis'); colorbar;

% ---- Label mask: bands exist for ALL times ----
maskTF = zeros(numel(t), numel(f), 'uint8');
idxNR   = abs(f - f0_NR)   <= BW_NR/2;
idxLTE  = abs(f - f0_LTE)  <= BW_LTE/2;
idxWLAN = abs(f - f0_WLAN) <= BW_WLAN/2;

maskTF(:, idxNR)   = 1;
maskTF(:, idxLTE)  = 2;
maskTF(:, idxWLAN) = 3;

% figure;
% imagesc(f/1e6, t*1e3, double(maskTF)); axis xy;
% xlim([-150 150]); ylim([0 Ttotal_ms]);
% xlabel('Frequency (MHz)'); ylabel('Time (ms)');
% title('Labeled spectrogram mask: 0=Noise, 1=NR, 2=LTE, 3=WLAN');
% colorbar;

%% ---- Time-domain plot of the composite signal xMix (amplitude vs time) ----
Nplot = numel(xMix);                       % plot full 0..Ttotal_ms duration
tSig  = (0:Nplot-1).'/FsCommon;            % time in seconds

% figure;
% plot(tSig*1e3, abs(xMix), 'LineWidth', 1.0);   % envelope magnitude
% grid on;
% xlabel('Time (ms)');
% ylabel('|xMix[n]|');
% title('Composite xMix magnitude (envelope) vs time');
% 
% % ( also plot I/Q (real/imag) for a short zoomed-in view
% zoom_ms = 0.2;                                % show first 0.2 ms
% Nz = min(Nplot, round(zoom_ms*1e-3*FsCommon));
% figure;
% plot((0:Nz-1).'/FsCommon*1e3, real(xMix(1:Nz)), 'LineWidth', 1.0); hold on;
% plot((0:Nz-1).'/FsCommon*1e3, imag(xMix(1:Nz)), 'LineWidth', 1.0);
% grid on;
% xlabel('Time (ms)');
% ylabel('Amplitude');
% legend('Re\{xMix\}','Im\{xMix\}', 'Location','best');
% title(sprintf('Composite xMix (I/Q) – first %.2f ms', zoom_ms));



%% ==========================================================
%  Robust occupied-band detection (forces 3 stripes robustly)
%  Inputs:
%    Sdb : (Nf x Nt) spectrogram in dB
%    f   : (Nf x 1) or (1 x Nf) frequency axis in Hz (monotone)
%  Outputs:
%    occBins : logical (Nf x 1) occupied bins after post-processing
%    bandsHz : [3 x 2] detected bands [fL fH] in Hz (top-3 by integrated power)
%
%  Key improvements vs your version:
%   - stronger gap filling (prevents one stripe splitting into two)
%   - mergeCloseBands() after bin->band conversion
%   - prune weak spurs by integrated band power
%   - keep exactly the top-3 strongest bands
%% ==========================================================

% --- 1) Time-averaged PSD per frequency bin (linear) ---
P = mean(10.^(Sdb/10), 2);        % Nf x 1
f = f(:);

% --- 2) Smooth across frequency (reduces speckle) ---
df = abs(f(2)-f(1));
smoothHz   = 2e6;                                  % smoothing span (Hz)
smoothBins = max(5, round(smoothHz/df));           % at least 5 bins
Psm = movmean(P, smoothBins);

% --- 3) Robust noise floor (median + MAD on low-power bins) ---
idxNoise = Psm <= prctile(Psm, 40);
muN  = median(Psm(idxNoise));
sigN = 1.4826 * mad(Psm(idxNoise), 1);             % robust scale (linear power)

% --- 4) Initial occupancy threshold ---
kOcc   = 6;                                        % higher => stricter (fewer spurs)
thrOcc = muN + kOcc*sigN;
occBins0 = (Psm >= thrOcc);

% --- 5) cleanup: fill gaps, remove small islands ---
minBandHz = 8e6;                                   % keep only reasonably wide bands
gapHz     = 8e6;                                   % IMPORTANT: fills OFDM notches so stripes don't split

minBandBins = max(1, round(minBandHz/df));
gapBins     = max(1, round(gapHz/df));

occBins = fillShortGapsBool(occBins0, gapBins);
occBins = removeShortRunsBool(occBins,  minBandBins);

% --- 6) Convert occupied bins -> candidate bands ---
candBandsHz = binsToBands(f, occBins);

% If nothing detected, exit
if isempty(candBandsHz)
    warning('No occupied bands detected. Consider lowering kOcc or minBandHz.');
    bandsHz = zeros(0,2);
    return;
end

% --- 7) Merge bands that are still too close (extra robustness) ---
mergeGapHz = 10e6;                                % merge close-by bands (Hz)
candBandsHz = mergeCloseBands(candBandsHz, mergeGapHz);

% --- 8) Rank by integrated band power and keep top-3 ---
bandPow = zeros(size(candBandsHz,1),1);
for k = 1:size(candBandsHz,1)
    idx = (f >= candBandsHz(k,1)) & (f <= candBandsHz(k,2));
    bandPow(k) = sum(Psm(idx));                   % integrated power over band
end

% Drop very weak bands (spurs) before taking top-3
minRelPow = 0.15;                                 % keep >=15% of strongest
keep = bandPow >= minRelPow*max(bandPow);
candBandsHz = candBandsHz(keep,:);
bandPow     = bandPow(keep);

% Now enforce exactly 3 bands (top-3 strongest)
Kwant = 3;
[~,ord] = sort(bandPow, 'descend');
ord = ord(1:min(Kwant, numel(ord)));
bandsHz = candBandsHz(ord,:);
bandsHz = sortrows(bandsHz, 1);

% Optional: rebuild occBins from final bands (so shading/printing matches)
occBins = false(size(f));
for k = 1:size(bandsHz,1)
    occBins = occBins | ((f >= bandsHz(k,1)) & (f <= bandsHz(k,2)));
end

% --- 9) Print exactly what was detected ---
fprintf('\n=== Final detected stripes/bands: %d ===\n', size(bandsHz,1));
for k = 1:size(bandsHz,1)
    fL = bandsHz(k,1); fH = bandsHz(k,2);
    fprintf('Stripe %d: fL=%+8.2f MHz, fH=%+8.2f MHz, BW=%6.2f MHz\n', ...
        k, fL/1e6, fH/1e6, (fH-fL)/1e6);
end
fprintf('=========================================\n\n');

% --- 10) Visualization (avg PSD + threshold + shaded final 3 bands) ---
figure;
plot(f/1e6, 10*log10(Psm+eps), 'LineWidth', 1.2); hold on;
yline(10*log10(thrOcc+eps), '--', 'LineWidth', 1.2);
grid on; xlabel('Frequency (MHz)'); ylabel('Avg power (dB)');
title('Robust occupied-band detection (keeps top-3 stripes)');

yl = ylim;
for k = 1:size(bandsHz,1)
    xpatch = [bandsHz(k,1) bandsHz(k,2) bandsHz(k,2) bandsHz(k,1)]/1e6;
    ypatch = [yl(1) yl(1) yl(2) yl(2)];
    patch(xpatch, ypatch, [0.8 0.8 0.8], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
end
uistack(findobj(gca,'Type','line'),'top');
legend('Avg PSD (smoothed)','Threshold','Top-3 bands','Location','best');







%% ==========================================================
%  From xMix + detected bandsHz -> isolate each stripe -> CP-corr scores
%  -> classify stripe as {NR, LTE, WLAN} based on highest score
%
%  Assumes we already have:
%   xMix (complex), FsCommon (Hz), bandsHz [Kx2] (Hz)
%  And  we still have from waveform generation:
%   FsWLAN, ofdmInfo, FsNR, FsLTE (or set them manually)
%% ==========================================================

% ---------- 0) Candidate lags (in samples at FsCommon) ----------
% WLAN Non-HT (20 MHz) uses NFFT=64
NfftWLAN = 64;
L_WLAN   = round(NfftWLAN * (FsCommon/FsWLAN));

% NR: use MATLAB's computed Nfft and sample rate
NfftNR = ofdmInfo.Nfft;
L_NR   = round(NfftNR * (FsCommon/FsNR));

% LTE: if LTE toolbox used, FsLTE is from infoLTE.SamplingRate.
% NfftLTE for 20 MHz LTE is typically 2048 (works for the fallback too).
NfftLTE = 1024;
L_LTE   = round(NfftLTE * (FsCommon/FsLTE));

lags = struct('NR',L_NR,'LTE',L_LTE,'WLAN',L_WLAN);
techList = {'NR','LTE','WLAN'};
lagList  = [lags.NR, lags.LTE, lags.WLAN];

fprintf('Candidate lags @ FsCommon=%.1f MHz: L_NR=%d, L_LTE=%d, L_WLAN=%d\n', ...
    FsCommon/1e6, lags.NR, lags.LTE, lags.WLAN);

% ---------- 1) CP-corr window parameters ----------
maxLag = max(lagList);
Nwin   = max(8*maxLag, round(0.5e-3*FsCommon));   % >= 8*maxLag or 0.5 ms
hop    = round(0.25*Nwin);                        % 25% hop
firOrd = 256;                                     % band-isolation filter order (tune 128..512)

% ---------- 2) Loop over detected stripes, isolate and classify ----------
K = size(bandsHz,1);
stripeLabels = strings(K,1);
stripeScores = zeros(K, numel(techList));

for k = 1:K
    band = bandsHz(k,:);  % [fL fH] in Hz

    % (A) isolate stripe in time series (mix down + LPF)
    [xBand, fcHz] = extractStripeTimeSeries(xMix, FsCommon, band(1), band(2), firOrd);

    % (B) CP-correlation score for each candidate lag (sliding, normalized)
    scores = zeros(1,numel(techList));
    for j = 1:numel(techList)
        L = lagList(j);
        [rhoVec, tRho] = slidingCpCorrCoeffFast(xBand, FsCommon, L, Nwin, hop);
        scores(j) = mean(rhoVec);     % summary score for this stripe vs this lag
    end

    stripeScores(k,:) = scores;

    % (C) pick technology with highest score
    [~,idxMax] = max(scores);
    stripeLabels(k) = techList{idxMax};

    fprintf('Stripe %d (%.1f..%.1f MHz, fc=%.1f MHz) ->  NR:%.3g  LTE:%.3g  WLAN:%.3g  => %s\n', ...
        k, band(1)/1e6, band(2)/1e6, fcHz/1e6, scores(1), scores(2), scores(3), stripeLabels(k));

    % (D) quick time-series plot of recovered stripe (envelope)
    figure;
    Nplot = min(numel(xBand), round(2e-3*FsCommon)); % first 2 ms
    tSig  = (0:Nplot-1).'/FsCommon;
    plot(tSig*1e3, abs(xBand(1:Nplot)), 'LineWidth', 1.0); grid on;
    xlabel('Time (ms)'); ylabel('|x_{band}[n]|');
    title(sprintf('Recovered Stripe %d (fc=%.1f MHz) envelope', k, fcHz/1e6));


    NfftSpec = 2^nextpow2(numel(xBand));      % or set fixed like 2^18 for speed
    X = fftshift(fft(xBand, NfftSpec));
    fSpec = (-NfftSpec/2:NfftSpec/2-1).' * (FsCommon/NfftSpec);  % baseband freq (Hz)
    
    Pdb = 10*log10(abs(X).^2 + 1e-12);
    Pdb = Pdb - max(Pdb);                    % normalize peak to 0 dB
    
    % Shift frequency axis back to the stripe center frequency
    fGlobal = fSpec + fcHz;                  % Hz
    
    figure;
    plot(fGlobal/1e6, Pdb, 'LineWidth', 1.1);
    grid on;
    xlabel('Frequency (MHz)');
    ylabel('Magnitude (dB, normalized)');
    title(sprintf('Stripe %d: FFT spectrum centered around f_c = %.1f MHz', k, fcHz/1e6));





    % (E)  plot normalized CP-corr curves for all tech lags
    figure; hold on; grid on;
    for j = 1:numel(techList)
        L = lagList(j);
        [rhoVec, tRho] = slidingCpCorrCoeffFast(xBand, FsCommon, L, Nwin, hop);
        plot(tRho*1e3, rhoVec, 'LineWidth', 1.0);
    end
    xlabel('Time (ms)'); ylabel('Normalized CP-corr coeff (windowed)');
    legend(techList, 'Location','best');
    title(sprintf('Stripe %d: CP-corr vs candidate lags', k));
end

% stripeLabels(k) now tells you which technology each stripe likely is.




