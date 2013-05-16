%***************************
%*  2D Radon transform synthetic example 
%*  
%*  2013,4,5
%*  References :
%*  
%*  
%***************************
%
%******************************
clear all; clc; close all;
datasetDirectory='../../dataset'
addpath src
addpath ../../ReadData3D_version1k/mha

t = 0:pi/20:2*pi;
R0 = 20; 
x0 = round(100/2); y0 = round(100/2);
xi = R0*cos(t)+x0;
yi = R0*sin(t)+y0;
roimask = poly2mask(xi,yi, 100,100);
pr_r = find(roimask);

subplot(3,3,1)
imshow(roimask)
 
theta = 0:180;
[R,xp] = radon(roimask,theta);
subplot(3,3,2)
imagesc(theta,xp,R); axis square;
title('R_{\theta} (X\prime)');
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
%colormap(hot);

subplot(3,3,3)
 [meanVHist,xout3]=hist(R,16);
 bar(xout3, meanVHist, 'b');
 grid on;

%% 
 t = 0:pi/20:2*pi;
R0 = 20; 
x0 = round(100/2); y0 = round(100/2);
xi = R0*cos(t)+x0+20;
yi = R0*sin(t)+y0;
roimask = poly2mask(xi,yi, 100,100);
pr_r = find(roimask);

subplot(3,3,4)
imshow(roimask)
 
theta = 0:180;
[R,xp] = radon(roimask,theta);
subplot(3,3,5)
imagesc(theta,xp,R); axis square;
title('R_{\theta} (X\prime)');
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
%colormap(hot);

subplot(3,3,6)
 [meanVHist,xout3]=hist(R,16);
 bar(xout3, meanVHist, 'b');
 grid on;
 
%% 
 t = 0:pi/20:2*pi;
R0 = 20; 
x0 = round(100/2); y0 = round(100/2);
xi = (R0+15)*cos(t)+x0+10;
yi = R0*sin(t)+y0;
roimask = poly2mask(xi,yi, 100,100);
pr_r = find(roimask);

subplot(3,3,7)
imshow(roimask)
 
theta = 0:180;
[R,xp] = radon(roimask,theta);
subplot(3,3,8)
imagesc(theta,xp,R); axis square;
title('R_{\theta} (X\prime)');
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
%colormap(hot);

subplot(3,3,9)
 [meanVHist,xout3]=hist(R,16);
 bar(xout3, meanVHist, 'b');
 grid on; 