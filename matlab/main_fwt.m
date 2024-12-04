clear;
load coeffs.mat;

img = double(imread('../images/airfield512x512.tif'));

% Compute all four filters by defination
LoD = db4;
n = 1:length(LoD);
HiD = -(-1).^(n-1) .* LoD(length(LoD) - n + 1);
LoR = LoD(length(LoD) - n + 1);
HiR = (-1).^(n-1) .* LoD(n);


% % Or equivalently
% [LoD, HiD, LoR, HiR] = wfilters('db4');

scales = 4;
[APPROXs, HORIZONTOLs, VERTICALs, DIAGONALs] = fwt(img, scales, LoD, HiD);

img_recon = ifwt(APPROXs, HORIZONTOLs, VERTICALs, DIAGONALs, scales, LoR, HiR);
figure
subplot(121)
imshow(uint8(img))
subplot(122)
imshow(uint8(img_recon))
mse(img, img_recon)
