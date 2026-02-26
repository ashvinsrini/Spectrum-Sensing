function maskRGB = maskToRGB_values_local(maskU8, labelVals, maskColors_u8)
% maskU8: HxW uint8 mask with values in labelVals (e.g., {0,63,127,191,6})
% labelVals: 1x5 uint8 [Noise NR LTE WLAN RADAR]
% maskColors_u8: 5x3 uint8 colors corresponding to labelVals order

[H,W] = size(maskU8);
maskRGB = zeros(H,W,3,'uint8');

for k = 1:numel(labelVals)
    v = labelVals(k);
    idx = (maskU8 == v);
    if any(idx,'all')
        for c = 1:3
            tmp = maskRGB(:,:,c);
            tmp(idx) = maskColors_u8(k,c);
            maskRGB(:,:,c) = tmp;
        end
    end
end
end
