%% Section 1 - Basis Functions
clc
clear all
close all
p = 65; %x
q = 65; %y
phase = atan(q/p);
disp('Phase Angle')
rad2deg(phase)
wavel =(2*pi)/(norm([p q]));
disp('Wavelength')
wavel
fftwave(p,q)
%% Section 2 - Linearity
clc
clear all
close all
F = [zeros(56,128); ones(16,128); zeros(56,128)];
G = F';
H = F + 2*G;
Fhat = fft2(F);
Ghat = fft2(G);
Hhat = fft2(H);
Fabsmax = max(abs(F(:)));
Gabsmax = max(abs(G(:)));
Habsmax = max(abs(H(:)));
HlogShift = log(1 + abs(fftshift(Hhat)));
figure(1)
subplot(3, 1, 1);
showgrey(F, 64, -Fabsmax, Fabsmax);
title('F')
subplot(3, 1, 2);
showgrey(G, 64, -Gabsmax, Gabsmax);
title('G')
subplot(3, 1, 3);
showgrey(H, 64, -Habsmax, Habsmax);
title('H')
figure(2)
subplot(3, 2, 1);
showgrey(log(1 + abs(fftshift(Fhat))));
title('F Spectra with log and fftshift')
subplot(3, 2, 3);
showgrey(log(1 + abs(fftshift(Ghat))));
title('G Spectra with log and fftshift')
subplot(3, 2, 5);
showgrey(log(1 + abs(fftshift(Hhat))));
title('H Spectra with log and fftshift')
subplot(3, 2, 2);
showgrey(fftshift(abs(Fhat)));
title('F Spectra without log')
subplot(3, 2, 4);
showgrey(fftshift(abs(Ghat)));
title('G Spectra without log')
subplot(3, 2, 6);
showgrey(fftshift(abs(Hhat)));
title('H Spectra without log')
figure(3)
subplot(2,1,1)
hist(Hhat)
title('Histogram before log')
subplot(2,1,2)
hist(HlogShift)
title('Histogram after log')
%% Multiplacation
clc
clear all
close all
F = [zeros(56,128); ones(16,128); zeros(56,128)];
sz = size(F,1);
G = F';
Fhat = fft2(F);
Ghat = fft2(G);
FG_fourier = fft2(F.*G);
FG_conv = (conv2(fftshift(Fhat),fftshift(Ghat),'same'))/(sz^2);
FG_conv = fftshift(FG_conv);
%figure(1)
%showgrey(F.*G);
figure(1)
subplot(2,1,1)
showfs(FG_fourier)
title('Multiplied in Spatial')
subplot(2,1,2)
showfs(FG_conv)
title('Convolved in Fourier')
%% Scaling
clc
clear all
close all
newf = [zeros(60, 128); ones(8, 128); zeros(60, 128)];
newg = [zeros(128, 48) ones(128, 32) zeros(128, 48)];
alpha = 45;
F = [newf .* ...
         newg];
Fhat = fft2(F);

figure(1)
subplot(1,2,1)
showgrey(newf)
title('New F')
subplot(1,2,2)
showgrey(newg)
title('New G')
figure(2)
subplot(1,2,1)
showgrey(F)
subplot(1,2,2)
showfs(Fhat)
G = rot(F, alpha);
Ghat = fft2(G);
Hhat = rot(fftshift(Ghat), -alpha);
figure(3)
subplot(1,3,1)
showgrey(G)
axis on
title(sprintf('%d degree Rotation in Spatial Domain', alpha))
subplot(1,3,2)
showfs(Ghat)
title(sprintf('%d degree Rotation in Frequency Domain', alpha))
axis on
subplot(1,3,3)
showgrey(log(1 + abs(Hhat)))
title(sprintf('%d degree Rotation in Frequency Domain', -alpha))
axis on
%% Rotations
clc
clear all
close all
newf = [zeros(60, 128); ones(8, 128); zeros(60, 128)];
newg = [zeros(128, 48) ones(128, 32) zeros(128, 48)];
alphas = [45 60 90];

F = [newf .* ...
         newg];
Fhat = fft2(F);
figure(1)
j = 0;
for i = 1:length(alphas)
alpha = alphas(i);
G = rot(F, alpha);
Ghat = fft2(G);
Hhat = rot(fftshift(Ghat), -alpha);
if(i == 1)
subplot(3,3,1)
showgrey(G)
axis on
title(sprintf('%d degree Rotation in Spatial Domain', alpha))
subplot(3,3,2)
showfs(Ghat)
title(sprintf('%d degree Rotation in Frequency Domain', alpha))
axis on
subplot(3,3,3)
showgrey(log(1 + abs(Hhat)))
title(sprintf('%d degree Rotation in Frequency Domain', -alpha))
axis on
elseif(i==2)
subplot(3,3,4)
showgrey(G)
axis on
title(sprintf('%d degree Rotation in Spatial Domain', alpha))
subplot(3,3,5)
showfs(Ghat)
title(sprintf('%d degree Rotation in Frequency Domain', alpha))
axis on
subplot(3,3,6)
showgrey(log(1 + abs(Hhat)))
title(sprintf('%d degree Rotation in Frequency Domain', -alpha))
axis on

else
subplot(3,3,7)
showgrey(G)
axis on
title(sprintf('%d degree Rotation in Spatial Domain', alpha))
subplot(3,3,8)
showfs(Ghat)
title(sprintf('%d degree Rotation in Frequency Domain', alpha))
axis on
subplot(3,3,9)
showgrey(log(1 + abs(Hhat)))
title(sprintf('%d degree Rotation in Frequency Domain', -alpha))
axis on
end
end
%% Pow2Im
clc
clear all
close all
a = 10^-9;
phone = phonecalc128;
few = few128;
nallo = nallo128;
pix = pow2image(phone,a);
pix1 = pow2image(few,a);
pix2 = pow2image(nallo,a);
img = randphaseimage(phone);
img1 = randphaseimage(few);
img2 = randphaseimage(nallo);
figure(1)
subplot(3,2,1)
showgrey((pix))
subplot(3,2,2)
showgrey(img)
subplot(3,2,3)
showgrey((pix1))
subplot(3,2,4)
showgrey(img1)
subplot(3,2,5)
showgrey((pix2))
subplot(3,2,6)
showgrey(img2)
%% Filtering procedure
% The Impulse response function
clc
clear all
close all
phone = phonecalc128;
n = 4;
var = [0.1 0.3 1 10 100];
figure = figure('Position',[100,100,1024,1200])
%figure(1)
j = 0;
for i= 1:length(var)
pix = gaussfft(deltafcn(128,128),var(i));
pix1 = discgaussfft(deltafcn(128,128),var(i));
subplot(5,2,i + j)
showgrey(pix)
title(sprintf('Gaussfft Variance %f', var(i)))
subplot(5,2,i+1+j)
showgrey(pix1)
title(sprintf('Discgaussfft Variance %f', var(i)))
j = i;
end
%% Smoothing
clc
close all
clear all
office = office256;
add = gaussnoise(office,16);
sap = sapnoise(office,0.1,255);
addg = gaussfft(add,0.5);
addmed = medfilt(add,5,5);
addlpf = ideal(add, 0.1);
addg1 = gaussfft(add,0.9);
addmed1 = medfilt(add,9,9);
addlpf1 = ideal(add, 0.01);
sapg = gaussfft(sap,0.5);
sapmed = medfilt(sap,5,5);
saplpf = ideal(sap, 0.1);
sapg1 = gaussfft(sap,0.9);
sapmed1 = medfilt(sap,9,9);
saplpf1 = ideal(sap, 0.01);
figure(1)
subplot(2,1,1)
showgrey(add)
subplot(2,1,2)
showgrey(sap)
figure(2)
subplot(3,3,1)
showgrey(add)
title('Original')
subplot(3,3,2)
showgrey(addg)
title('Gaussian Smoothing Variance 0.5')
subplot(3,3,3)
showgrey(addg1)
title('Gaussian Smoothing Variance 0.9')
subplot(3,3,4)
showgrey(add)
title('Original')
subplot(3,3,5)
showgrey(addmed)
title('Median Filter Window 5')
subplot(3,3,6)
showgrey(addmed1)
title('Median Filter Window 9')
subplot(3,3,7)
showgrey(add)
title('Original')
subplot(3,3,8)
showgrey(addlpf)
title('Low pass filter cuttoff 0.1')
subplot(3,3,9)
showgrey(addlpf1)
title('Low pass filter cuttoff 0.01')
figure(3)
subplot(3,3,1)
showgrey(sap)
title('Original')
subplot(3,3,2)
showgrey(sapg)
title('Gaussian Smoothing Variance 0.5')
subplot(3,3,3)
showgrey(sapg1)
title('Gaussian Smoothing Variance 0.9')
subplot(3,3,4)
showgrey(sap)
title('Original')
subplot(3,3,5)
showgrey(sapmed)
title('Median Filter Window 5')
subplot(3,3,6)
showgrey(sapmed1)
title('Median Filter Window 9')
subplot(3,3,7)
showgrey(sap)
title('Original')
subplot(3,3,8)
showgrey(saplpf)
title('Low pass filter cuttoff 0.1')
subplot(3,3,9)
showgrey(saplpf1)
title('Low pass filter cuttoff 0.01')
%% Subsampling
clc
close all
clear all
img = phonecalc256;
smoothimg = img;
N=5;
figure(1)
%title('Gaussian Smoothing Variance 0.3')
for i=1:N
if i>1    % generate subsampled versions
   img = rawsubsample(img);
   smoothimg = ideal(smoothimg, 0.1);
   smoothimg = rawsubsample(smoothimg);
end
subplot(2, N, i)
showgrey(img)
subplot(2, N, i+N)
showgrey(smoothimg)
end