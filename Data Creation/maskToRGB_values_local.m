
function rgb = maskToRGB_values_local(maskU8, labelVals, colors_u8)
% Map mask values (e.g., {0,63,127,191,6}) to RGB colors for visualization.
% labelVals and colors_u8 must be in order [Noise NR LTE WLAN RADAR].

H = size(maskU8,1); W = size(maskU8,2);
rgb = zeros(H,W,3,'uint8');

% default -> Noise index
idx = ones(H,W,'uint8');

for k = 1:numel(labelVals)
    idx(maskU8 == labelVals(k)) = uint8(k);
end

rgb(:,:,1) = reshape(colors_u8(double(idx),1), H, W);
rgb(:,:,2) = reshape(colors_u8(double(idx),2), H, W);
rgb(:,:,3) = reshape(colors_u8(double(idx),3), H, W);
end