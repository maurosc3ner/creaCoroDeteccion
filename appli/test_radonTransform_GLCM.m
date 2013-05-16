%***************************
%*  2D Radon transform example
%*  
%*  2013,4,5
%*  References :
%*  
%*  
%***************************
%
% 
% 
% 
%
% 
%******************************
clear all; clc; close all;
datasetDirectory='../../dataset'
addpath src
addpath ../../ReadData3D_version1k/mha

addpath ../src/GLCM_Features4

occOff=round(4);
occuDir=[0 occOff;   
        -occOff occOff;    
        -occOff 0;    
        -occOff -occOff];    
%         0 -occOff;     
%         occOff -occOff;    
%         occOff 0;     
%         occOff occOff]
meanContr=[];
meanHomog=[];
meanEnerg=[];
meanEntro=[];

%% file lecture
cprFilename=strcat(datasetDirectory,'/','dt07_z106_S.mhd');

info = mha_read_header(cprFilename);
V = mha_read_volume(info);
[dims]=size(V);
[minmaxV]=[min(min(min(V))) max(max(max(V)))];

%%bounding box 4mm
radio=4;
R0 = radio/info.PixelDimensions(1); 

I=V(round(dims(1)/2)-R0:round(dims(1)/2)+R0,...
        round(dims(2)/2)-R0:round(dims(2)/2)+R0);

%I=V(round(size(V,1)/4):round(3*size(V,1)/4),...
%        round(size(V,1)/4):round(3*size(V,1)/4));
%I = zeros(100,100);
%I(25:75, 25:75) = 1;
figure;
subplot(2,2,1)
imagesc(I),axis image;

theta = 0:180;
[R,xp] = radon(I,theta);
subplot(2,2,2)
imagesc(theta,xp,R);
title('R_{\theta} (X\prime)');
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
%colormap(hot);
colorbar
[mean(R(:)) std(R(:))]


%% co-occurrence
[GLCM1 SI] = graycomatrix(R,'Offset',occuDir,'NumLevels',32);
stats = graycoprops(GLCM1,{'contrast','homogeneity','energy',});
     stats.Contrast;
     stats.Homogeneity;
     stats.Energy;
contr=mean(stats.Contrast);
homog=mean(stats.Homogeneity);
energ=mean(stats.Energy);
entro=entropy(R);
meanContr=[meanContr contr];
meanHomog=[meanHomog homog];
meanEntro=[meanEntro entro];
meanEnerg=[meanEnerg energ];

hIntProf=sum(R,1);
vIntProf=sum(R,2);
vIntProf=vIntProf./100000;
vIntProf=vIntProf+abs(min(vIntProf));

figure;
subplot(2,1,1)
plot(hIntProf)
subplot(2,1,2)
plot(vIntProf)
std(vIntProf)
mean(vIntProf)
hold all
plot(round(77/2),vIntProf(round(77/2)),'*g')
%pr_r=find(vIntProf>mean(vIntProf))
pr_r2=find(1:77>round(77/2-2*std(vIntProf)) & 1:77<round(77/2+2*std(vIntProf)))
plot(pr_r2,vIntProf(pr_r2),'*r')
hold off

t = 0:pi/20:2*pi;
R0 = 20; 
x0 = round(100/2); y0 = round(100/2);
xi = R0*cos(t)+x0;
yi = R0*sin(t)+y0;
roimask = poly2mask(xi,yi, 100,100);
pr_r = find(roimask);

subplot(2,2,3)
imshow(roimask)
% [meanVHist,xout3]=hist(R,32);
% bar(xout3, meanVHist, 'b'); 
[R,xp] = radon(roimask,theta);
subplot(2,2,4)
imagesc(theta,xp,R);
title('R_{\theta} (X\prime)');
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
%colormap(hot);
colorbar

[mean(R(:)) std(R(:))];

minmax=[min(min(R)) max(max(R))]
%[GLCM1] = graycomatrix(R,'Offset',occuDir,'NumLevels',64,'Symmetric',true);
[GLCM1 SI] = graycomatrix(R,'GrayLimits',minmax,'Offset',occuDir,'NumLevels',64,'Symmetric',true);
    
stats = graycoprops(GLCM1,{'contrast','homogeneity','energy',});
     stats.Contrast
     stats.Homogeneity
     stats.Energy
 contr=mean(stats.Contrast);
 homog=mean(stats.Homogeneity);
 energ=mean(stats.Energy);
 entro=entropy(R);

stats4=GLCM_Features4(GLCM1,1)

meanContr=[meanContr contr];
meanHomog=[meanHomog homog];
meanEntro=[meanEntro entro];
meanEnerg=[meanEnerg energ];

hIntProf=sum(R,1);
vIntProf=sum(R,2);
figure;
subplot(2,1,1)
plot(hIntProf)
subplot(2,1,2)
plot(vIntProf)
std(vIntProf);