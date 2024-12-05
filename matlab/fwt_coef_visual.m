function fwt_coef_visual(scales, img_height, img_width, APPROXs, HORIZONTOLs, VERTICALs, DIAGONALs, margin_enable)
    %% Find maximal and minimal values
    min_val = Inf;
    max_val = -Inf;

    for scale_to_plot = 1:scales
        approx = APPROXs{scale_to_plot};
        horizontal = HORIZONTOLs{scale_to_plot};
        vertical = VERTICALs{scale_to_plot};
        diagonal = DIAGONALs{scale_to_plot};

        % Update max and min values
        max_val = max(max_val, max(approx(:)));
        max_val = max(max_val, max(horizontal(:))); 
        max_val = max(max_val, max(vertical(:)));
        max_val = max(max_val, max(diagonal(:)));

        min_val = min(min_val, min(approx(:)));
        min_val = min(min_val, min(horizontal(:))); 
        min_val = min(min_val, min(vertical(:)));
        min_val = min(min_val, min(diagonal(:)));
    end
    
    %% Plot FWT coefficients spliced
    spiral_plot = fwt_coef_splice(scales, img_height, img_width, APPROXs, HORIZONTOLs, VERTICALs, DIAGONALs, max_val, margin_enable);
    figure;
    disp(size(spiral_plot));
    imagesc(spiral_plot, [min_val, max_val]); % Keep brightness range fixed
    colormap gray;
    % colorbar
    title(sprintf('Wavelet Coefficients (with white margin)'));
    axis square;
    axis off;

    %% Plot each scale
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
        imagesc(spiral_plot, [min_val, max_val]); % Keep brightness range fixed
        colormap gray;
        % colorbar
        title(sprintf('Wavelet Coefficients (Scale %d)', scale_to_plot));
        axis square;
        axis off;
    end
end
