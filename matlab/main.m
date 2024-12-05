clear;
close all;

% mode = 'low bit rates';
mode = 'all bit rates';

[dct_bit_rates, dct_psnr_vec] = main_dct(mode);
[fwt_bit_rates, fwt_psnr_vec] = main_fwt(mode);

%%
figure;
plot(dct_bit_rates, dct_psnr_vec, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(fwt_bit_rates, fwt_psnr_vec, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
grid minor;
set(gca, 'FontSize', 10, 'GridAlpha', 0.3); 
title('PSNR vs Bit Rate', 'FontSize', 14);
xlabel('Bit Rate (bits per pixel)', 'FontSize', 12);
ylabel('PSNR (dB)', 'FontSize', 12);
legend("DCT","FWT", 'Location', 'best', 'FontSize', 10); 
set(gca, 'FontSize', 10, 'GridAlpha', 0.3); 
hold off;