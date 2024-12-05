function [bit_rates, psnr_vec] = main_dct(mode)
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
    step_size = 1; % you can change this value
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
    if mode == 'all bit rates'
        step_range = round(2.^(0:9));
    elseif mode == 'low bit rates'
        step_range = round(2.^(6:0.3:9)); % Step size range
    else
        error('Please enter a valid mode')
    end
    bit_rates = []; % Stores the average bit rate across all 8x8 DCT coefficients for different step sizes
    psnr_vec = []; % Stores the average PSNR values for different step sizes
    % 
    % step_sample = 2.^[0 5 9];
    % boats_sample = {boats};

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
        psnr_step = PSNR([boats, harbour, peppers], [boats_recon, harbour_recon, peppers_recon]); % PSNR for boats image
    %     harbour_psnr = PSNR(harbour, harbour_recon); % PSNR for harbour image
    %     peppers_psnr = PSNR(peppers, peppers_recon); % PSNR for peppers image
    %     psnr_step = (boats_psnr + harbour_psnr + peppers_psnr) / 3; % Average PSNR for this step size
        psnr_vec = [psnr_vec, psnr_step]; % Store average PSNR for this step size

        % Save Samples for plotting
        % if ismember(step, step_sample)
        %     boats_sample = [boats_sample boats_recon];
        % end
    end

    % figure;
    % for i = 1:numel(boats_sample)
    %     subplot(1, numel(boats_sample), i)
    %     imshow(uint8(boats_sample{i}))
    %     if i == 1
    %         title('Original')
    %     else
    %         title(sprintf('Reconstructed, step = %d', step_sample(i-1)))
    %     end
    % end


%     figure;
%     plot(bit_rates, psnr_vec, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
%     grid on;
%     title('PSNR vs Bit Rate (DCT)', 'FontSize', 14);
%     xlabel('Bit Rate (bits per pixel)', 'FontSize', 12);
%     ylabel('PSNR (dB)', 'FontSize', 12);
%     legend('PSNR Curve', 'Location', 'best', 'FontSize', 10); 
%     set(gca, 'FontSize', 10, 'GridAlpha', 0.3); 

end

function [bit_rates, psnr_vec] = main_fwt(mode)
% clear
    load coeffs.mat;

    %% Compute all four filters by defination
    LoD = db4;
    n = 1:length(LoD);
    HiD = -(-1).^(n-1) .* LoD(length(LoD) - n + 1);
    LoR = LoD(length(LoD) - n + 1);
    HiR = (-1).^(n-1) .* LoD(n);
    % % Or equivalently
    % [LoD, HiD, LoR, HiR] = wfilters('db4');


    files = ["../images/harbour512x512.tif"; "../images/boats512x512.tif"; "../images/peppers512x512.tif"];

    %% rate-PSNR curve, Sec 3.3
    if mode == 'all bit rates'
        step_range = round(2.^(0:9));
    elseif mode == 'low bit rates'
        step_range = round(2.^(6:0.3:9)); % Step size range
    else
        error('Please enter a valid mode')
    end
    bit_rates = []; % Stores the average bit rate across all 8x8 DCT coefficients for different step sizes
    psnr_vec = []; % Stores the average PSNR values for different step sizes
    img_mse= zeros(10, 3);
    coeff_mse = zeros(10, 3);
    
    % boats_sample = {double(imread(files(2)))};
    % step_sample = 2.^[0 5 9];

    scales = 4;
    for steps = step_range
        SUBBANDS = {};
        if steps ==1 || steps ==16 || steps == 512
            fig=figure();
        end
        imgs=[];
        img_recons=[];
        for img_idx=1:3 % loop through files
            %% 3.2 FWT
            img =double(imread(files(img_idx)));

            [APPROXs, HORIZONTOLs, VERTICALs, DIAGONALs] = fwt(img, scales, LoD, HiD);

            if img_idx==1 && steps ==1 % We only visualize FWT coef for the 1st image. 
                % Combine coefficients in a spiral layout for visualization
                [img_height, img_width]=size(img);
                fwt_coef_visual(scales, img_height, img_width, APPROXs, HORIZONTOLs, VERTICALs, DIAGONALs, 1)
            end

            %% 3.3 Quantization
            APPROXs_quant = quantizerCell(APPROXs, steps);
            HORIZONTOLs_quant = quantizerCell(HORIZONTOLs, steps);
            VERTICALs_quant = quantizerCell(VERTICALs, steps);
            DIAGONALs_quant = quantizerCell(DIAGONALs, steps);

            %% 3.4 IFWT Reconstruction, Distortion, and BitRate
            img_recon = ifwt(APPROXs_quant, HORIZONTOLs_quant, VERTICALs_quant, DIAGONALs_quant, scales, LoR, HiR);

            %% Plot original&reconstructed image
            if steps ==1 || steps ==16 || steps == 512 % We only Plot original&reconstructed image for 3 steps 
                ax1=subplot(3, 2, img_idx*2-1,'Parent', fig);
                imshow(uint8(img), 'Parent', ax1);
                str=sprintf("Original Image, Step = %d", steps);
                title(ax1, str);
                ax2=subplot(3, 2, img_idx*2, 'Parent', fig);
                imshow(uint8(img_recon), 'Parent', ax2);
                str=sprintf("Reconstructed Image, Step = %d", steps);
                title(ax2, str);
            end

            % if img_idx == 2                              
            %     if ismember(steps, step_sample)
            %         boats_sample = [boats_sample img_recon];
            %     end
            % end
            
            %% Compare mse between DWT coefficients and images
            mse_1 = 0;  mse_2 = 0;  mse_3 = 0;  mse_4 = 0;
            j=1;
            deno = size(APPROXs{scales} ,1) * size(APPROXs{scales} ,2);
            mse_1 = mse_1 + mse(APPROXs{scales}, APPROXs_quant{scales}) * size(APPROXs{scales} ,1) * size(APPROXs{scales} ,2);
            SUBBAND={};
            SUBBAND{1}=APPROXs_quant{scales};
            for i = 1:numel(APPROXs)
                j=j+1;
                SUBBAND{j}=HORIZONTOLs_quant{i};
                mse_2 = mse_2 + mse(HORIZONTOLs{i}, HORIZONTOLs_quant{i}) * size(APPROXs{i} ,1) * size(APPROXs{i} ,2) ;

                j=j+1;
                SUBBAND{j}=VERTICALs_quant{i};
                mse_3 = mse_3 + mse(VERTICALs{i}, VERTICALs_quant{i}) * size(APPROXs{i} ,1) * size(APPROXs{i} ,2) ;

                j=j+1;
                SUBBAND{j}=DIAGONALs_quant{i};
                mse_4 = mse_4 + mse(DIAGONALs{i}, DIAGONALs_quant{i}) * size(APPROXs{i} ,1) * size(APPROXs{i} ,2) ;
                deno = deno + 3* size(APPROXs{i} ,1) * size(APPROXs{i} ,2);
            end
            img_mse(steps, img_idx) = mse(img, img_recon);
            coeff_mse(steps, img_idx) = (mse_1 + mse_2 + mse_3 + mse_4) / deno;
            fprintf("step = %d, img = %s: \t img_mse = %f, coeff_mse = %f\n", steps, files(img_idx), img_mse(steps, img_idx), coeff_mse(steps, img_idx))
            imgs=[imgs, img];
            img_recons=[img_recons, img_recon];
            if isempty(SUBBANDS)
                SUBBANDS = SUBBAND;
            else
                for i=1:j
                    SUBBANDS{i}= [SUBBANDS{i},SUBBAND{i}];
                end
            end
            %spiral_plot = fwt_coef_splice(scales, img_height, img_width, APPROXs, HORIZONTOLs, VERTICALs, DIAGONALs,max_val, 0);

        end
        psnr_step = PSNR(imgs,img_recons);
        psnr_vec=[psnr_vec,psnr_step];

        bit_rate_step = 0;
        subband_total_size=0;
        for i=1:numel(SUBBANDS)
            bit_rate_step = bit_rate_step + computeBitRate(SUBBANDS{i}) * size(SUBBANDS{i}, 1)*size(SUBBANDS{i}, 2);
            subband_total_size=subband_total_size + size(SUBBANDS{i}, 1)*size(SUBBANDS{i}, 2);
        end
        bit_rate_step = bit_rate_step/subband_total_size;

        bit_rates =[bit_rates, bit_rate_step]; % Store average bit rate for this step size
    end
%     figure;
%     plot(bit_rates, psnr_vec, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
%     grid on;
%     title('PSNR vs Bit Rate (FWT)', 'FontSize', 14);
%     xlabel('Bit Rate (bits per pixel)', 'FontSize', 12);
%     ylabel('PSNR (dB)', 'FontSize', 12);
%     legend('PSNR Curve', 'Location', 'best', 'FontSize', 10); 
%     set(gca, 'FontSize', 10, 'GridAlpha', 0.3); 

    % figure;
    % for i = 1:numel(boats_sample)
    %     subplot(1, numel(boats_sample), i)
    %     imshow(uint8(boats_sample{i}))
    %     if i == 1
    %         title('Original')
    %     else
    %         title(sprintf('Reconstructed, step = %d', step_sample(i-1)))
    %     end
    % end
end