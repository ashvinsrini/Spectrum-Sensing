function rgb8 = SdbToRGB_local(Sdb, cmap, useFixedCLim, clim_dB)
% Convert Sdb (Nf x Nt) -> RGB image like imagesc(..., Sdb.')
Splot = double(Sdb.');  % Nt x Nf

if useFixedCLim
    vmin = clim_dB(1);
    vmax = clim_dB(2);
else
    vmax = max(Splot(:));
    vmin = vmax - 80;
end

Sclip = min(max(Splot, vmin), vmax);
idx0  = round(255 * (Sclip - vmin) / max(vmax - vmin, eps)); % 0..255
idx   = uint8(idx0 + 1);                                     % 1..256
rgb   = ind2rgb(idx, cmap);                                  % double [0,1]
rgb8  = im2uint8(rgb);
end