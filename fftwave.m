function fftwave(u, v, sz)
if (nargin < 2)
  error('Requires at least two input arguments.')
end
if (nargin == 2)
sz = 128; end
Fhat = zeros(sz);
Fhat(u, v) = 1;
F = ifft2(Fhat);
Fabsmax = max(abs(F(:)));
subplot(3, 2, 1);
imshow(Fhat);
title(sprintf('Fhat: (u, v) = (%d, %d)', u, v))
% What is done by these instructions?
if (u <= sz/2)
uc = u - 1; else
  uc = u - 1 - sz;
end
if (v <= sz/2)
  vc = v - 1;
else
  vc = v - 1 - sz;
end
wavelength = 0.0;  % Replace by correct expression
amplitude = 0.0;   % Replace by correct expression
subplot(3, 2, 2);
imshow(fftshift(Fhat));
title(sprintf('centered Fhat: (uc, vc) = (%d, %d)', uc, vc))
subplot(3, 2, 3);
imshow(real(F),  [-Fabsmax Fabsmax]);
title('real(F)')
subplot(3, 2, 4);
imshow(imag(F), [-Fabsmax Fabsmax]);
title('imag(F)')
subplot(3, 2, 5);
imshow(abs(F), [-Fabsmax Fabsmax]);
title(sprintf('abs(F) (amplitude %f)', amplitude))
subplot(3, 2, 6);
imshow(angle(F),  [-pi pi]);
title(sprintf('angle(F) (wavelength %f)', wavelength))