% Sec 2.1
% Function to compute the DCT transform matrix A of size MxM
function matA = matA_DCT(M)
    matA = zeros(M, M); % Initialize the matrix
    for i = 0:M-1
        for k = 0:M-1
            if i == 0
                a = sqrt(1 / M); % Compute alpha_0
            else
                a = sqrt(2 / M); % Compute alpha_i for i > 0
            end
            % Populate the matrix element A[i, k]
            matA(i+1, k+1) = a * cos( (2*k+1)*i*pi / (2*M) );
        end
    end
end
