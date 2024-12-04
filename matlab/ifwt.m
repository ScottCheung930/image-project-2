function reconstructed_img = ifwt(APPROXs, HORIZONTOLs, VERTICALs, DIAGONALs, scales, LoR, HiR)
    % To recover the image, only the smallest approxiamation is needed
    APPROX_tmp = APPROXs{end};
    for i = scales : -1 : 1
        APPROX_tmp = filterbank_synthesis_2d(APPROX_tmp, HORIZONTOLs{i}, VERTICALs{i}, DIAGONALs{i}, LoR, HiR);
    end
    reconstructed_img = APPROX_tmp;
end

function reconstruction_2d = filterbank_synthesis_2d(APPROX, HORIZONTOL, VERTICAL, DIAGONAL, LoR, HiR)
    % W ~ width of the current coeff; H ~ height...
    [W, H] = size(APPROX);

    % Filter the vertical axis
    item_Lo = zeros(H*2, W);
    item_Hi = zeros(H*2, W);
    for i = 1:W
        item_Lo(:, i) = filterbank_synthesis_1d(APPROX(:, i)', HORIZONTOL(:, i)', LoR, HiR);
        item_Hi(:, i) = filterbank_synthesis_1d(VERTICAL(:, i)', DIAGONAL(:, i)', LoR, HiR);
    end

    % Filter the horizontal axis
    reconstruction_2d = zeros(H*2, W*2);
    for i = 1:2*H
        reconstruction_2d(i, :) = filterbank_synthesis_1d(item_Lo(i, :), item_Hi(i, :), LoR, HiR);
    end

end

function reconstruction_1d = filterbank_synthesis_1d(low_band, high_band, LoR, HiR)
    L = length(LoR);
    
    % Up Sampling
    low_band = upsample(low_band, 2);
    high_band = upsample(high_band, 2);

    % Add circular padding to avoid marginal effect
    pad_low_band = padarray(low_band, [0, L-1], "circular");
    pad_high_band = padarray(high_band, [0, L-1], "circular");

    % Filtering using two filter banks
    pad_filtered_Lo = conv(pad_low_band, LoR, "full");
    pad_filtered_Hi = conv(pad_high_band, HiR, "full");
    
    % Remove padding, the resulted data should have the same length with
    % the original data
    filtered_Lo = pad_filtered_Lo( L+1 : end-(L-1)*2+1 );
    filtered_Hi = pad_filtered_Hi( L+1 : end-(L-1)*2+1 );

    reconstruction_1d = filtered_Lo + filtered_Hi;
end

