function [bit_rates, psnr_vec] = main_fwt()
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

    step_range = 2.^(0:9); % Step size range
    bit_rates = []; % Stores the average bit rate across all 8x8 DCT coefficients for different step sizes
    psnr_vec = []; % Stores the average PSNR values for different step sizes
    img_mse= zeros(10, 3);
    coeff_mse = zeros(10, 3);

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

            %% Compare mse between DWT coefficients and images
            mse_1 = 0;  mse_2 = 0;  mse_3 = 0;  mse_4 = 0;
            j=1;
            deno = size(APPROXs{scales} ,1) * size(APPROXs{scales} ,2);
            mse_1 = mse_1 + mse(APPROXs{scales}, APPROXs_quant{scales}) * size(APPROXs{scales} ,1) * size(APPROXs{scales} ,2);
            SUBBAND={};
            SUBBAND{1}=APPROXs{scales};
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
        for i=1:(2^scales - (scales -1))
            bit_rate_step = bit_rate_step + computeBitRate(SUBBANDS{i}) * size(SUBBANDS{i}, 1)*size(SUBBANDS{i}, 2);
            subband_total_size=subband_total_size+size(SUBBANDS{i}, 1)*size(SUBBANDS{i}, 2);
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
end