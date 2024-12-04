clear;
load coeffs.mat;
img = double(imread('../images/harbour512x512.tif'));

%% Compute all four filters by defination
LoD = db4;
n = 1:length(LoD);
HiD = -(-1).^(n-1) .* LoD(length(LoD) - n + 1);
LoR = LoD(length(LoD) - n + 1);
HiR = (-1).^(n-1) .* LoD(n);
% % Or equivalently
% [LoD, HiD, LoR, HiR] = wfilters('db4');

%% FWT
scales = 4;
[APPROXs, HORIZONTOLs, VERTICALs, DIAGONALs] = fwt(img, scales, LoD, HiD);

%% Combine coefficients in a spiral layout for visualization
for scale_to_plot = 1:scales
    approx = APPROXs{scale_to_plot};
    horizontal = HORIZONTOLs{scale_to_plot};
    vertical = VERTICALs{scale_to_plot};
    diagonal = DIAGONALs{scale_to_plot};
    % Get the size of the coefficients
    [m, n] = size(approx);    
    % Create a larger matrix to accommodate the spiral layout
    spiral_plot = zeros(2*m, 2*n);    
    % Place the coefficients in the respective quadrants
    spiral_plot(1:m, 1:n) = approx;                    % Top-left: Approximation
    spiral_plot(1:m, n+1:2*n) = horizontal;            % Top-right: Horizontal
    spiral_plot(m+1:2*m, 1:n) = vertical;              % Bottom-left: Vertical
    spiral_plot(m+1:2*m, n+1:2*n) = diagonal;          % Bottom-right: Diagonal    
    figure;
    imagesc(spiral_plot);
    colormap gray(256);
    % colorbar;
    title(sprintf('Wavelet Coefficients (Scale %d)', scale_to_plot));
    axis square;
    axis off;
end

%% Quantization
steps = 1;
APPROXs_quant = quantizerCell(APPROXs, steps);
HORIZONTOLs_quant = quantizerCell(HORIZONTOLs, steps);
VERTICALs_quant = quantizerCell(VERTICALs, steps);
DIAGONALs_quant = quantizerCell(DIAGONALs, steps);

%% IFWT
img_recon = ifwt(APPROXs_quant, HORIZONTOLs_quant, VERTICALs_quant, DIAGONALs_quant, scales, LoR, HiR);

%% Plot original&reconstructed image
figure
subplot(121)
imshow(uint8(img))
title('Original Image')
subplot(122)
imshow(uint8(img_recon))
title('Reconstructed Image')

%% Compare mse between DWT coefficients and images
mse_1 = 0;
mse_2 = 0;
mse_3 = 0;
mse_4 = 0;
deno = 0;
for i = 1:numel(APPROXs)
    mse_1 = mse_1 + mse(APPROXs{i}, APPROXs_quant{i}) * size(APPROXs{i} ,1) * size(APPROXs{i} ,2);
    mse_2 = mse_2 + mse(HORIZONTOLs{i}, HORIZONTOLs_quant{i}) * size(APPROXs{i} ,1) * size(APPROXs{i} ,2) ;
    mse_3 = mse_3 + mse(VERTICALs{i}, VERTICALs_quant{i}) * size(APPROXs{i} ,1) * size(APPROXs{i} ,2) ;
    mse_4 = mse_4 + mse(DIAGONALs{i}, DIAGONALs_quant{i}) * size(APPROXs{i} ,1) * size(APPROXs{i} ,2) ;
    deno = deno + size(APPROXs{i} ,1) * size(APPROXs{i} ,2);
end
img_mse = mse(img, img_recon)
coeff_mse = (mse_1 + mse_2 + mse_3 + mse_4) / 4 / deno
