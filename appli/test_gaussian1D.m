%% Example of gaussian equation and filter
clc;close all; clear all;

sigma = 10;
size = 6*sigma+1;

%% 1D
x = linspace(-size / 2, size / 2, size);        %range
gaussFilter = (1/(sqrt(2*pi)*sigma))*exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normal
LoG=(1-(x.^2)./sigma^2).*exp(-x .^ 2 / (2 * sigma ^ 2));

figure;
subplot(2,1,1)
plot(gaussFilter),title('Gaussian with sigma =10')
subplot(2,1,2)
plot(LoG),title('LoG with sigma =10')

%procesar algo por un gaussiano
y = rand(500,1);
plot(1:500,y)
yfilt = filter (gaussFilter,1, y);
plot(1:500,yfilt)
hold on;
yfilt = conv (y, gaussFilter, 'same');
plot(1:500,yfilt)
y = -10:50;
yfilt = conv (y, gaussFilter, 'same');
plot(-10:50,yfilt)

%% 2D
[X Y] =meshgrid(-64:64,-64:64);
gauss2d=(1/(2*pi*sigma^2))*exp(-(X.^2+Y.^2)/(2*sigma^2));
gauss2d=gauss2d./sum(gauss2d(:));
LoG2d=(2-(X.^2+Y.^2)./sigma^2).*exp(-(X.^2+Y.^2)/(2*sigma^2));

figure;
subplot(1,2,1)
surf(gauss2d),title('Gaussian with sigma =10')
subplot(1,2,2)
surf(LoG2d),title('LoG with sigma =10')

[fX fY]=gradient(gauss2d);
levelset=1./(1+abs(fX)+abs(fY));

[fX fY]=gradient(LoG2d);
levelset=1./(1+abs(fX)+abs(fY));

