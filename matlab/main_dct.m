clear;

%% Print out the values of the matrix A, Sec 2.1
M = 8;
matA = matA_DCT(8);
disp(matA)

%% Plot a graph of the quantizer function, Sec 2.2
original_vec = linspace(-10, 10, 1000);
quant_vec = midTreadQuant(original_vec, 1);
figure;
plot(original_vec, '--', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);
hold on;
plot(quant_vec, '-', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
title('Comparison of Original and Quantized Data', 'FontSize', 14);
xlabel('Index', 'FontSize', 12);
ylabel('Amplitude', 'FontSize', 12);
legend('Original Data', 'Quantized Data', 'Location', 'best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 10, 'GridAlpha', 0.3);
axis tight;
set(gcf, 'Color', 'w');

% Load images
boats = double(imread("../images/boats512x512.tif"));
harbour = double(imread("../images/harbour512x512.tif"));
peppers = double(imread("../images/peppers512x512.tif"));
[lenX, lenY] = size(boats);

%% Compare d with the mean squared error between the original and the quantized DCT coefficients, Sec 2.3
step_size = 10; % you can change this value
img = harbour; % you can choose another image to test
img_dct_mat = []; % matrix that stores the DCT coeffs
img_dct_quant_mat = []; % matrix that stores the quantized DCT coeffs
img_recon = zeros(lenX, lenY);
for i = 1 : lenX/M
    for k = 1 : lenY/M
        img_block = img(M*(i-1)+1 : M*i, M*(k-1)+1 : M*k);
        img_block_dct = dct2(img_block); % DCT
        img_dct_mat = [img_dct_mat; img_block_dct];
        img_block_dct_quant = midTreadQuant(img_block_dct, step_size); % quantize the DCT coeffs
        img_dct_quant_mat = [img_dct_quant_mat; img_block_dct_quant];
        img_recon(M*(i-1)+1 : M*i, M*(k-1)+1 : M*k) = idct2(img_block_dct_quant); % IDCT
    end
end
d = mse(img, img_recon); % MSE of original/reconstructed images
d_dct = mse(img_dct_mat, img_dct_quant_mat); % MSE of original/quantized DCT coeffs
disp("The difference between MSE of DCT Coefficients and MSE of Images = ")
disp(mae(d-d_dct))

%% rate-PSNR curve, Sec 2.3
step_range = 2.^(0:9); % Step size range
bit_rates = []; % Stores the average bit rate across all 8x8 DCT coefficients for different step sizes
psnr_vec = []; % Stores the average PSNR values for different step sizes

% Loop over different quantizer step sizes
for step = step_range
    block_dct_quant_mat = []; % Stores quantized DCT coefficients for all positions in all blocks
    boats_recon = zeros(lenX, lenY); % Reconstructed image for boats
    harbour_recon = zeros(lenX, lenY); % Reconstructed image for harbour
    peppers_recon = zeros(lenX, lenY); % Reconstructed image for peppers

    % Loop through all 8x8 blocks of the images
    for i = 1 : lenX/M
        for k = 1 : lenY/M
            % Extract and quantize DCT coefficients for each block
            boats_block_dct_quant = midTreadQuant( dct2(boats(M*(i-1)+1 : M*i, M*(k-1)+1 : M*k)), step );
            harbour_block_dct_quant = midTreadQuant( dct2(harbour(M*(i-1)+1 : M*i, M*(k-1)+1 : M*k)), step );
            peppers_block_dct_quant = midTreadQuant( dct2(peppers(M*(i-1)+1 : M*i, M*(k-1)+1 : M*k)), step );

            % Flatten and store the quantized coefficients for all images
            block_dct_quant = [boats_block_dct_quant(:), harbour_block_dct_quant(:), peppers_block_dct_quant(:)];
            block_dct_quant_mat = [block_dct_quant_mat, block_dct_quant];

            % Perform inverse DCT to reconstruct the block
            boats_recon(M*(i-1)+1 : M*i, M*(k-1)+1 : M*k) = idct2(boats_block_dct_quant);
            harbour_recon(M*(i-1)+1 : M*i, M*(k-1)+1 : M*k) = idct2(harbour_block_dct_quant);
            peppers_recon(M*(i-1)+1 : M*i, M*(k-1)+1 : M*k) = idct2(peppers_block_dct_quant);
        end
    end

    % Calculate the bit rate for each coefficient position across all blocks
    N = size(block_dct_quant_mat, 1); % Number of DCT coefficients (64 in this case)
    bit_rate_N = []; % Stores bit rates for all coefficient positions
    for n = 1:N
        bit_rate_n = computeBitRate(block_dct_quant_mat(n, :)); % Compute the bit rate for position n
        bit_rate_N = [bit_rate_N, bit_rate_n]; % Store bit rate for this position
    end
    bit_rate_step = mean(bit_rate_N); % Average bit rate across all coefficient positions
    bit_rates = [bit_rates, bit_rate_step]; % Store average bit rate for this step size

    % Calculate PSNR for each image and take the average
    boats_psnr = PSNR(boats, boats_recon); % PSNR for boats image
    harbour_psnr = PSNR(harbour, harbour_recon); % PSNR for harbour image
    peppers_psnr = PSNR(peppers, peppers_recon); % PSNR for peppers image
    psnr_step = (boats_psnr + harbour_psnr + peppers_psnr) / 3; % Average PSNR for this step size
    psnr_vec = [psnr_vec, psnr_step]; % Store average PSNR for this step size
end


figure;
plot(bit_rates, psnr_vec, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
title('PSNR vs Bit Rate', 'FontSize', 14);
xlabel('Bit Rate (bits per pixel)', 'FontSize', 12);
ylabel('PSNR (dB)', 'FontSize', 12);
legend('PSNR Curve', 'Location', 'best', 'FontSize', 10); 
set(gca, 'FontSize', 10, 'GridAlpha', 0.3); 

