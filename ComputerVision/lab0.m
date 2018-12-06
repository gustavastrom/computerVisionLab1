%%
close all
clear all
clc
nallo = nallo256float;
negtransf = (255: -1: 0)';
neg3 = 255 - nallo;
% neg3 = compose(negtransf,round(nallo)+1);
% figure(1)
% showgrey(nallo);
% figure(2)
% showgrey(log(0.05 + nallo));
subplot(3,2,1)
showgrey(nallo)
subplot(3,2,2)
hist(nallo(:))
subplot(3,2,3)
showgrey(neg3)
subplot(3,2,4)
hist(neg3(:))
subplot(3,2,5)
showgrey(log(0.0 + nallo))
subplot(3,2,6)
showgrey(log(0.1 + nallo))
%%
clc
close all
clear all
Canoe = 0;
load ('canoe256.mat', 'Canoe')
phone = phonecalc256;

negtransf = (255: -1: 0)';
neg3 = compose(negtransf,Canoe+1);
%showgrey(log(phone))
neg1 = -Canoe;
plus0 = Canoe;
plus1 = Canoe + 1;
figure(1)
%showgrey(neg1)
hist(neg1(:))
neg2 = 255 - Canoe;
diff = neg3 - neg2;
figure(2)
hist(neg2(:))
%showgrey(plus0)
figure(3)
%showgrey(plus1)
hist(neg3(:))
%% Color tables
clc
clear all
close all
Canoe = 0;
load ('canoe256.mat', 'Canoe')
figure(1)
image(Canoe+1)
negcolormapcol = linspace(1,0,256)';
colormap([negcolormapcol negcolormapcol negcolormapcol])
figure(2)
showgrey(Canoe,linspace(1,0,256),0,255)