close all
clc
clear
Fhat = zeros(128,128);
p = 5;
q = 9;
fftwave(p,q,128)
% Fhat(p,q) = 1;
% F = ifft2(Fhat);
% Fabsmax = max(abs(F(:)));
% figure
% imshow(real(F),  [-Fabsmax Fabsmax])
% figure
% imshow(imag(F), [-Fabsmax Fabsmax])
% figure
% imshow(abs(F), [-Fabsmax Fabsmax])
% figure
% imshow(angle(F), [-pi pi])
%imshow(Fhat);