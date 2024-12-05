function spiral_plot = fwt_coef_splice(scales, img_height, img_width, APPROXs, HORIZONTOLs, VERTICALs, DIAGONALs, max_val, margin_enable)
    if margin_enable == 1
        spiral_plot = zeros(img_height+2^(scales+1), img_width+2^(scales+1));
    else
        spiral_plot = zeros(img_height, img_width);
    end
    
    for scale_to_plot = 1:scales
        approx = APPROXs{scale_to_plot};
        horizontal = HORIZONTOLs{scale_to_plot};
        vertical = VERTICALs{scale_to_plot};
        diagonal = DIAGONALs{scale_to_plot};
        % Get the size of the coefficients
        [m, n] = size(approx);
        if margin_enable ==1
            margin = 2^(scales-scale_to_plot);
        else
            margin =0;
        end
        % Place the coefficients in the respective quadrants
        if scale_to_plot == scales
            spiral_plot(1:m+2*margin, 1: n+2*margin)=max_val.*ones(m+2*margin,  n+2*margin);
            spiral_plot(1+margin:m+margin, 1+margin:n+margin) = approx;                    % Top-left: Approximation
        end
        spiral_plot(1:(m+2*margin), (n+2*margin+1): (2*(n+2*margin)))=max_val.*ones(m+2*margin,  n+2*margin);
        spiral_plot(1+margin:m+margin, n+1+3*margin:2*n+3*margin) = horizontal;            % Top-right: Horizontal

        spiral_plot(m+2*margin+1:2*(m+2*margin), 1:n+2*margin)=max_val.*ones(m+2*margin,  n+2*margin);
        spiral_plot(m+1+3*margin:2*m+3*margin, 1+margin:n+margin) = vertical;              % Bottom-left: Vertical

        spiral_plot(m+2*margin+1:2*(m+2*margin), n+2*margin+1: 2*(n+2*margin))=max_val.*ones(m+2*margin,  n+2*margin);
        spiral_plot(m+1+3*margin:2*m+3*margin, n+1+3*margin:2*n+3*margin) = diagonal;          % Bottom-right: Diagonal
        %fprintf("Scale = %d, [m, n] = [%d, %d] ", scale_to_plot, m,n);
        %fprintf("\n");
    end
    
end