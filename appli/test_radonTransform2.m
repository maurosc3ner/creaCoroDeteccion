clear all; clc; close all;
datasetDirectory='../../dataset'
addpath src
addpath ../src/GLCM_Features4
addpath ../../ReadData3D_version1k/mha


meanRadon=[];
stdRadon=[];

meanEntr=[];

wsplot=4;
hsplot=4;
windex=[5,6,1,2,9,10,13,14,3,4,7,8,11,12,15,16];
theta = 0:180;
t = 0:pi/20:2*pi;
R0 = 20; 
x0 = round(100/2); y0 = round(100/2);
xi = R0*cos(t)+x0;
yi = R0*sin(t)+y0;
roimask = poly2mask(xi,yi, 100,100);
pr_r = find(roimask);

%% HEALTHY

%% file lecture
cprFilename=strcat(datasetDirectory,'/','dt07_z106_S.mhd');

info = mha_read_header(cprFilename)
V = mha_read_volume(info);
[dims]=size(V);

%%bounding box 4mm
radio=4;
R0 = radio/info.PixelDimensions(1); 

I=V(round(dims(1)/2)-R0:round(dims(1)/2)+R0,...
        round(dims(2)/2)-R0:round(dims(2)/2)+R0);
newMax = 1;
newMin = 0;
    
figure;
subplot(hsplot,wsplot,windex(1))
imagesc(I),axis image;title('dt07 106 S.mhd');

[counts, x]=hist(I(:),16);
bar(x, counts,'r')
meanEntr=[meanEntr -sum(counts.*log2(counts))];


[R,xp] = radon(I,theta);
meanRadon=[meanRadon mean(R(:))];
stdRadon=[stdRadon std(R(:))];
subplot(hsplot,wsplot,windex(2))
imagesc(theta,xp,R);axis square;
title(strcat(num2str(mean(R(:))),'|',num2str(std(R(:)))));
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
colormap(gray);

subplot(hsplot,wsplot,windex(3))
imshow(roimask)
% [meanVHist,xout3]=hist(R,16);
% bar(xout3, meanVHist, 'b'); 
[R,xp] = radon(roimask,theta);
subplot(hsplot,wsplot,windex(4))
imagesc(theta,xp,R);axis square;
title('R_{\theta} (X\prime)');
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);

%% file lecture
cprFilename=strcat(datasetDirectory,'/','dt00_z40.mhd');

info = mha_read_header(cprFilename);
V = mha_read_volume(info);
[dims]=size(V);

%%bounding box 4mm
radio=4;
R0 = radio/info.PixelDimensions(1); 

I=V(round(dims(1)/2)-R0:round(dims(1)/2)+R0,...
        round(dims(2)/2)-R0:round(dims(2)/2)+R0);

subplot(hsplot,wsplot,windex(5))
imagesc(I),axis image;title('dt00 40 S.mhd');

[counts, x]=hist(I(:),16);
bar(x, counts,'r')
meanEntr=[meanEntr -sum(counts.*log2(counts))];

[R,xp] = radon(I,theta);
meanRadon=[meanRadon mean(R(:))];
stdRadon=[stdRadon std(R(:))];
subplot(hsplot,wsplot,windex(6))
imagesc(theta,xp,R);axis square;
title(strcat(num2str(mean(R(:))),'|',num2str(std(R(:)))));
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
colormap(gray);


%% file lecture
cprFilename=strcat(datasetDirectory,'/','dt08_z105_S.mhd');

info = mha_read_header(cprFilename);
V = mha_read_volume(info);
[dims]=size(V);

%%bounding box 4mm
radio=4;
R0 = radio/info.PixelDimensions(1); 

I=V(round(dims(1)/2)-R0:round(dims(1)/2)+R0,...
        round(dims(2)/2)-R0:round(dims(2)/2)+R0);


subplot(hsplot,wsplot,windex(7))
imagesc(I),axis image;title('dt08 105 S.mhd');

[counts, x]=hist(I(:),16);
bar(x, counts,'r')
meanEntr=[meanEntr -sum(counts.*log2(counts))];

[R,xp] = radon(I,theta);
meanRadon=[meanRadon mean(R(:))];
stdRadon=[stdRadon std(R(:))];
subplot(hsplot,wsplot,windex(8))
imagesc(theta,xp,R);axis square;
title(strcat(num2str(mean(R(:))),'|',num2str(std(R(:)))));
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
colormap(gray);

%% MIXED  

%% file lecture
cprFilename=strcat(datasetDirectory,'/','dt08_z68_M.mhd');

info = mha_read_header(cprFilename);
V = mha_read_volume(info);
[dims]=size(V);

%%bounding box 4mm
radio=4;
R0 = radio/info.PixelDimensions(1); 

I=V(round(dims(1)/2)-R0:round(dims(1)/2)+R0,...
        round(dims(2)/2)-R0:round(dims(2)/2)+R0);
    
subplot(hsplot,wsplot,windex(9))
imagesc(I),axis image;title('dt08 68 M.mhd');

[counts, x]=hist(I(:),16);
bar(x, counts,'r')
meanEntr=[meanEntr -sum(counts.*log2(counts))];

[R,xp] = radon(I,theta);
meanRadon=[meanRadon mean(R(:))];
stdRadon=[stdRadon std(R(:))];
subplot(hsplot,wsplot,windex(10))
imagesc(theta,xp,R);axis square;
title(strcat(num2str(mean(R(:))),'|',num2str(std(R(:)))));
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
colormap(gray);


%% CALCIFIED
%% file lecture
cprFilename=strcat(datasetDirectory,'/','dt07_z78_C.mhd');

info = mha_read_header(cprFilename);
V = mha_read_volume(info);
[dims]=size(V);

%%bounding box 4mm
radio=4;
R0 = radio/info.PixelDimensions(1); 

I=V(round(dims(1)/2)-R0:round(dims(1)/2)+R0,...
        round(dims(2)/2)-R0:round(dims(2)/2)+R0);

subplot(hsplot,wsplot,windex(11))
imagesc(I),axis image;title('dt07 78 C.mhd');

[counts, x]=hist(I(:),16);
bar(x, counts,'r')
meanEntr=[meanEntr -sum(counts.*log2(counts))]

[R,xp] = radon(I,theta);
meanRadon=[meanRadon mean(R(:))];
stdRadon=[stdRadon std(R(:))];
subplot(hsplot,wsplot,windex(12))
imagesc(theta,xp,R);axis square;
title(strcat(num2str(mean(R(:))),'|',num2str(std(R(:)))));
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
colormap(gray);


%% file lecture
cprFilename=strcat(datasetDirectory,'/','dt08_z32_C.mhd');

info = mha_read_header(cprFilename);
V = mha_read_volume(info);
[dims]=size(V);

%%bounding box 4mm
radio=4;
R0 = radio/info.PixelDimensions(1); 

I=V(round(dims(1)/2)-R0:round(dims(1)/2)+R0,...
        round(dims(2)/2)-R0:round(dims(2)/2)+R0);

subplot(hsplot,wsplot,windex(13))
imagesc(I),axis image;title('dt08 32 C.mhd');

[counts, x]=hist(I(:),16);
bar(x, counts,'r')
meanEntr=[meanEntr -sum(counts.*log2(counts))]

[R,xp] = radon(I,theta);
meanRadon=[meanRadon mean(R(:))];
stdRadon=[stdRadon std(R(:))];
subplot(hsplot,wsplot,windex(14))
imagesc(theta,xp,R);axis square;
title(strcat(num2str(mean(R(:))),'|',num2str(std(R(:)))));
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
colormap(gray);


%% file lecture
cprFilename=strcat(datasetDirectory,'/','dt08_z51_C.mhd');

info = mha_read_header(cprFilename);
V = mha_read_volume(info);
[dims]=size(V);

%%bounding box 4mm
radio=4;
R0 = radio/info.PixelDimensions(1); 

I=V(round(dims(1)/2)-R0:round(dims(1)/2)+R0,...
        round(dims(2)/2)-R0:round(dims(2)/2)+R0);

subplot(hsplot,wsplot,windex(15))
imagesc(I),axis image;title('dt08 51 C.mhd');

[counts, x]=hist(I(:),16);
bar(x, counts,'r')
meanEntr=[meanEntr -sum(counts.*log2(counts))];

[R,xp] = radon(I,theta);
meanRadon=[meanRadon mean(R(:))];
stdRadon=[stdRadon std(R(:))];
subplot(hsplot,wsplot,windex(16))
imagesc(theta,xp,R);axis square;
title(strcat(num2str(mean(R(:))),'|',num2str(std(R(:)))));
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
colormap(gray);

%% Feature Curves
figure;
curvesWS=2;
indC=5:7;
indM=4;
indS=1:3;
%% mean radon transform
subplot(curvesWS,curvesWS-1,1)
plot(0:6,meanRadon),title('mean RT');
hold all;
sP=plot(indS-1,meanRadon(indS),'*c');
set(sP, 'Color',[0 1 0]);
cP=plot(indC-1,meanRadon(indC),'*r');
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indM-1,meanRadon(indM),'*m');
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Mean');
hleg1 = legend('Mean','Healthy','Calcified','Mixed');
%hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;

%% std radon transform
subplot(curvesWS,curvesWS-1,2)
plot(0:6,stdRadon),title('std RT');
hold all;
sP=plot(indS-1,stdRadon(indS),'*c');
set(sP, 'Color',[0 1 0]);
cP=plot(indC-1,stdRadon(indC),'*r');
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indM-1,stdRadon(indM),'*m');
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','std');
hleg1 = legend('std','Healthy','Calcified','Mixed');
%hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;
