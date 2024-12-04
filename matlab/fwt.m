function [APPROXs, HORIZONTOLs, VERTICALs, DIAGONALs] = fwt(img, scales, LoD, HiD)
    APPROXs = {};
    HORIZONTOLs = {};
    VERTICALs = {};
    DIAGONALs = {};

    for i = 1:scales
        [APPROX, HORIZONTOL, VERTICAL, DIAGONAL] = filterbank_analysis_2d(img, LoD, HiD);
        APPROXs = [APPROXs APPROX];
        HORIZONTOLs = [HORIZONTOLs HORIZONTOL];
        VERTICALs = [VERTICALs VERTICAL];
        DIAGONALs = [DIAGONALs DIAGONAL];

        % Only retain the approximations in each iteration
        img = APPROX;
    end
end

function [APPROX, HORIZONTOL, VERTICAL, DIAGONAL] = filterbank_analysis_2d(data_2d, LoD, HiD)
    % W ~ width of original img; H ~ height...
    [W, H] = size(data_2d);
    
    % Filter the horizontal axis
    intm_Lo = zeros(H, floor(W/2));
    intm_Hi = zeros(H, floor(W/2));
    for i = 1:H
        [intm_Lo(i, :), intm_Hi(i, :)] = filterbank_analysis_1d(data_2d(i, :), LoD, HiD);
    end

    % Filter the vertical axis
    APPROX = zeros(floor(H/2), floor(W/2)); % low-low
    HORIZONTOL = zeros(floor(H/2), floor(W/2)); % low-high
    VERTICAL = zeros(floor(H/2), floor(W/2)); % high-low
    DIAGONAL = zeros(floor(H/2), floor(W/2)); % high-high
    for i = 1:floor(W/2)
        [APPROX(:, i), HORIZONTOL(:, i)] = filterbank_analysis_1d(intm_Lo(:, i)', LoD, HiD);
        [VERTICAL(:, i), DIAGONAL(:, i)] = filterbank_analysis_1d(intm_Hi(:, i)', LoD, HiD);
    end
end

function [low_band, high_band] = filterbank_analysis_1d(data_1d, LoD, HiD)
    L = length(LoD);

    % Add circular padding to avoid marginal effect
    pad_data_1d = padarray(data_1d, [0, L-1], "circular");

    % Filtering using two filter banks
    pad_filtered_Lo = conv(pad_data_1d, LoD, "full");
    pad_filtered_Hi = conv(pad_data_1d, HiD, "full");

    % Remove padding, the resulted data should have the same length with
    % the original data
    filtered_Lo = pad_filtered_Lo( (L-1)*2 : (L-1)*2+length(data_1d)-1 );
    filtered_Hi = pad_filtered_Hi( (L-1)*2 : (L-1)*2+length(data_1d)-1 );

    % Down Sampling
    low_band = downsample(filtered_Lo, 2);
    high_band = downsample(filtered_Hi, 2);
end