function fwt_coef_visual(scales, img_height, img_width, APPROXs, HORIZONTOLs, VERTICALs, DIAGONALs, margin_enable)
    
    
    
    %% find maximal value
    for scale_to_plot = 1:scales
        approx = APPROXs{scale_to_plot};
        horizontal = HORIZONTOLs{scale_to_plot};
        vertical = VERTICALs{scale_to_plot};
        diagonal = DIAGONALs{scale_to_plot};
        max_val=max(max(approx));
        max_val = max(max_val, max(max(horizontal))); 
        max_val = max(max_val, max(max(vertical)));
        max_val = max(max_val, max(max(diagonal)));
    end
    
    %% plot FWT coeffecient spliced
    
    spiral_plot = fwt_coef_splice(scales, img_height, img_width, APPROXs, HORIZONTOLs, VERTICALs, DIAGONALs, max_val, margin_enable);
    figure;
    disp(size(spiral_plot))
    imagesc(spiral_plot);
    colormap gray(256);
    % colorbar;
    title(sprintf('Wavelet Coefficients (with white margin)'));
    axis square;
    axis off;
    %% plot each scale
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

end