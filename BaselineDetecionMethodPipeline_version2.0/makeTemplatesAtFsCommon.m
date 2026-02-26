function tmpl = makeTemplatesAtFsCommon(FsCommon, Ntotal, carrierNR)
% Precompute baseband templates @ FsCommon length Ntotal (single)
    n = (0:Ntotal-1).';

    % ---- WLAN template (Non-HT CBW20) ----
    cfgWLAN = wlanNonHTConfig('ChannelBandwidth','CBW20','MCS',4,'PSDULength',1000);
    bitsWiFi = randi([0 1], cfgWLAN.PSDULength*8, 1);
    txWLAN   = wlanWaveformGenerator(bitsWiFi, cfgWLAN);
    FsWLAN   = wlanSampleRate(cfgWLAN);
    xWLAN = resampleAny(txWLAN, FsWLAN, FsCommon);
    xWLAN = repToLen(xWLAN, Ntotal);

    % ---- NR template (random full grid) ----
    ofdmInfo = nrOFDMInfo(carrierNR);
    FsNR = ofdmInfo.SampleRate;
    Nsc  = carrierNR.NSizeGrid * 12;
    Nsym = carrierNR.SymbolsPerSlot;
    M = 16;
    qamSymbols = qammod(randi([0 M-1], Nsc*Nsym, 1), M, 'UnitAveragePower', true);
    gridNR = reshape(qamSymbols, [Nsc Nsym 1]);
    txNR = nrOFDMModulate(carrierNR, gridNR);
    xNR = resampleAny(txNR(:), FsNR, FsCommon);
    xNR = repToLen(xNR, Ntotal);

    % ---- LTE-like OFDM template (simple, CP-style) ----
    % (Works even without LTE toolbox; sufficient for baseline CP-lag scoring)
    FsLTE = 30.72e6;
    Nfft  = 2048;
    NscLTE = 1200;       % 20 MHz LTE active subcarriers
    NsymLTE = 14;        % one subframe worth, then repeat
    xLTE0 = simpleOFDMLike(FsLTE, Nfft, NscLTE, NsymLTE);
    xLTE  = resampleAny(xLTE0(:), FsLTE, FsCommon);
    xLTE  = repToLen(xLTE, Ntotal);

    % ---- RADAR pulse train template (baseband LFM pulse train) ----
    xRAD = radarPulseTrainBaseband(FsCommon, Ntotal);

    tmpl.NR    = single(xNR(:));
    tmpl.LTE   = single(xLTE(:));
    tmpl.WLAN  = single(xWLAN(:));
    tmpl.RADAR = single(xRAD(:));
end