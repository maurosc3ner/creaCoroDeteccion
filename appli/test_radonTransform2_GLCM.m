clear all; clc; close all;
datasetDirectory='../../dataset'
addpath src
addpath ../src/GLCM_Features4

occOff=round(2);

% %una direccion
% occuDir=[0 occOff;   
%          0 2*occOff;
%          0 3*occOff;
%          0 4*occOff]; 

%dos direcciones (45,0)
% occuDir=[0 occOff;   
%         -occOff occOff];

% 8 direcciones    
occuDir=[0 occOff;
        -occOff occOff;
        -occOff 0;
        -occOff -occOff];


nLevels=64;
meanAutoc=[];
meanCorrp=[];
meanContr=[];
meanHomog=[];
meanEnerg=[];
meanEntrop=[];
meanVar=[];
meanCS=[];
meanCP=[];
meanMP=[];
meanRadon=[];
stdRadon=[];

wsplot=3;
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
    
%I = (I - minmaxV(1))*(newMax - newMin)/(minmaxV(2) - minmaxV(1)) + newMin;

% figure;imagesc(I),axis image,colormap gray;
% figure;hist(I);
% figure;imagesc(Inorm),axis image,colormap gray;
% figure;hist(Inorm);
%I=V(round(size(V,1)/4):round(3*size(V,1)/4),...
%        round(size(V,1)/4):round(3*size(V,1)/4));
%I = zeros(100,100);
%I(25:75, 25:75) = 1;
figure;
subplot(2,wsplot,1)
imagesc(I),axis image;title('dt07 106 S.mhd');

theta = 0:180;
[R,xp] = radon(I,theta);
[mean(R(:)) std(R(:))];
subplot(2,wsplot,2)
imagesc(theta,xp,R);axis square;
title(strcat(num2str(mean(R(:))),'|',num2str(std(R(:)))));
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
colormap(gray);

[minmaxV]=[min(min(min(V))) max(max(max(V)))];
[minmaxR]=[min(R(:)) max(R(:))];
%% co-occurrence
[GLCM1 SI] = graycomatrix(I,'GrayLimits',minmaxV,'Offset',occuDir,'NumLevels',nLevels,'Symmetric',true);
stats4=GLCM_Features4(GLCM1,1);
meanAutoc=[meanAutoc mean(stats4.autoc)];
meanCorrp=[meanCorrp mean(stats4.corrp)];
meanContr=[meanContr mean(stats4.contr)];
meanHomog=[meanHomog mean(stats4.homop)];
meanEnerg=[meanEnerg mean(stats4.energ)];
meanEntrop=[meanEntrop mean(stats4.entro)];
meanVar=[meanVar mean(stats4.sosvh)];
meanCS=[meanCS mean(stats4.cshad)];
meanCP=[meanCP mean(stats4.cprom)];
meanMP=[meanMP mean(stats4.maxpr)];
meanRadon=[meanRadon mean(R(:))];
stdRadon=[stdRadon std(R(:))];

% hIntProf=sum(R,1);
% vIntProf=sum(R,2);
% subplot(2,wsplot,3)
% plot(hIntProf)
% subplot(2,wsplot,6);
% plot(vIntProf);
% std(vIntProf);


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
%I = (I - minmaxV(1))*(newMax - newMin)/(minmaxV(2) - minmaxV(1)) + newMin;

%I=V(round(size(V,1)/4):round(3*size(V,1)/4),...
%        round(size(V,1)/4):round(3*size(V,1)/4));
%I = zeros(100,100);
%I(25:75, 25:75) = 1;
figure;
subplot(2,wsplot,1)
imagesc(I),axis image;title('dt00 40 S.mhd');

theta = 0:180;
[R,xp] = radon(I,theta);
[mean(R(:)) std(R(:))];
subplot(2,wsplot,2)
imagesc(theta,xp,R);axis square;
title(strcat(num2str(mean(R(:))),'|',num2str(std(R(:)))));
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
%colormap(hot);

[minmaxV]=[min(min(min(V))) max(max(max(V)))];
[minmaxR]=[min(R(:)) max(R(:))];
%% co-occurrence
[GLCM1 SI] = graycomatrix(I,'GrayLimits',minmaxV,'Offset',occuDir,'NumLevels',nLevels,'Symmetric',true);
stats4=GLCM_Features4(GLCM1,1);
meanAutoc=[meanAutoc mean(stats4.autoc)];
meanCorrp=[meanCorrp mean(stats4.corrp)];
meanContr=[meanContr mean(stats4.contr)];
meanHomog=[meanHomog mean(stats4.homop)];
meanEnerg=[meanEnerg mean(stats4.energ)];
meanEntrop=[meanEntrop mean(stats4.entro)];
meanVar=[meanVar mean(stats4.sosvh)];
meanCS=[meanCS mean(stats4.cshad)];
meanCP=[meanCP mean(stats4.cprom)];
meanMP=[meanMP mean(stats4.maxpr)];
meanRadon=[meanRadon mean(R(:))];
stdRadon=[stdRadon std(R(:))];

% hIntProf=sum(R,1);
% vIntProf=sum(R,2);
% subplot(2,wsplot,3)
% plot(hIntProf)
% subplot(2,wsplot,6);
% plot(vIntProf);
% std(vIntProf);


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
%I = (I - minmaxV(1))*(newMax - newMin)/(minmaxV(2) - minmaxV(1)) + newMin;

%I=V(round(size(V,1)/4):round(3*size(V,1)/4),...
%        round(size(V,1)/4):round(3*size(V,1)/4));
%I = zeros(100,100);
%I(25:75, 25:75) = 1;
figure;
subplot(2,wsplot,1)
imagesc(I),axis image;title('dt08 105 S.mhd');

theta = 0:180;
[R,xp] = radon(I,theta);
[mean(R(:)) std(R(:))];
subplot(2,wsplot,2)
imagesc(theta,xp,R);axis square;
title(strcat(num2str(mean(R(:))),'|',num2str(std(R(:)))));
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
%colormap(hot);

[minmaxV]=[min(min(min(V))) max(max(max(V)))];
[minmaxR]=[min(R(:)) max(R(:))];
%% co-occurrence
[GLCM1 SI] = graycomatrix(I,'GrayLimits',minmaxV,'Offset',occuDir,'NumLevels',nLevels,'Symmetric',true);
stats4=GLCM_Features4(GLCM1,1);
meanAutoc=[meanAutoc mean(stats4.autoc)];
meanCorrp=[meanCorrp mean(stats4.corrp)];
meanContr=[meanContr mean(stats4.contr)];
meanHomog=[meanHomog mean(stats4.homop)];
meanEnerg=[meanEnerg mean(stats4.energ)];
meanEntrop=[meanEntrop mean(stats4.entro)];
meanVar=[meanVar mean(stats4.sosvh)];
meanCS=[meanCS mean(stats4.cshad)];
meanCP=[meanCP mean(stats4.cprom)];
meanMP=[meanMP mean(stats4.maxpr)];
meanRadon=[meanRadon mean(R(:))];
stdRadon=[stdRadon std(R(:))];
% hIntProf=sum(R,1);
% vIntProf=sum(R,2);
% subplot(2,wsplot,3)
% plot(hIntProf)
% subplot(2,wsplot,6);
% plot(vIntProf);
% std(vIntProf);

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
%I = (I - minmaxV(1))*(newMax - newMin)/(minmaxV(2) - minmaxV(1)) + newMin;

%I=V(round(size(V,1)/4):round(3*size(V,1)/4),...
%        round(size(V,1)/4):round(3*size(V,1)/4));
%I = zeros(100,100);
%I(25:75, 25:75) = 1;
figure;
subplot(2,wsplot,1)
imagesc(I),axis image;title('dt08 68 M.mhd');

theta = 0:180;
[R,xp] = radon(I,theta);
[mean(R(:)) std(R(:))];
subplot(2,wsplot,2)
imagesc(theta,xp,R);axis square;
title(strcat(num2str(mean(R(:))),'|',num2str(std(R(:)))));
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
%colormap(hot);


[minmaxV]=[min(min(min(V))) max(max(max(V)))];
[minmaxR]=[min(R(:)) max(R(:))];
%% co-occurrence
[GLCM1 SI] = graycomatrix(I,'GrayLimits',minmaxV,'Offset',occuDir,'NumLevels',nLevels,'Symmetric',true);
stats4=GLCM_Features4(GLCM1,1);
meanAutoc=[meanAutoc mean(stats4.autoc)];
meanCorrp=[meanCorrp mean(stats4.corrp)];
meanContr=[meanContr mean(stats4.contr)];
meanHomog=[meanHomog mean(stats4.homop)];
meanEnerg=[meanEnerg mean(stats4.energ)];
meanEntrop=[meanEntrop mean(stats4.entro)];
meanVar=[meanVar mean(stats4.sosvh)];
meanCS=[meanCS mean(stats4.cshad)];
meanCP=[meanCP mean(stats4.cprom)];
meanMP=[meanMP mean(stats4.maxpr)];
meanRadon=[meanRadon mean(R(:))];
stdRadon=[stdRadon std(R(:))];
% hIntProf=sum(R,1);
% vIntProf=sum(R,2);
% subplot(2,wsplot,3)
% plot(hIntProf)
% subplot(2,wsplot,6);
% plot(vIntProf);
% std(vIntProf);


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
%%I = (I - minmaxV(1))*(newMax - newMin)/(minmaxV(2) - minmaxV(1)) + newMin;

%I=V(round(size(V,1)/4):round(3*size(V,1)/4),...
%        round(size(V,1)/4):round(3*size(V,1)/4));
%I = zeros(100,100);
%I(25:75, 25:75) = 1;
figure;
subplot(2,wsplot,1)
imagesc(I),axis image;title('dt07 78 C.mhd');

theta = 0:180;
[R,xp] = radon(I,theta);
[mean(R(:)) std(R(:))];
subplot(2,wsplot,2)
imagesc(theta,xp,R);axis square;
title(strcat(num2str(mean(R(:))),'|',num2str(std(R(:)))));
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
%colormap(hot);

[minmaxV]=[min(min(min(V))) max(max(max(V)))];
[minmaxR]=[min(R(:)) max(R(:))];
%% co-occurrence
[GLCM1 SI] = graycomatrix(I,'GrayLimits',minmaxV,'Offset',occuDir,'NumLevels',nLevels,'Symmetric',true);
stats4=GLCM_Features4(GLCM1,1);
meanAutoc=[meanAutoc mean(stats4.autoc)];
meanCorrp=[meanCorrp mean(stats4.corrp)];
meanContr=[meanContr mean(stats4.contr)];
meanHomog=[meanHomog mean(stats4.homop)];
meanEnerg=[meanEnerg mean(stats4.energ)];
meanEntrop=[meanEntrop mean(stats4.entro)];
meanVar=[meanVar mean(stats4.sosvh)];
meanCS=[meanCS mean(stats4.cshad)];
meanCP=[meanCP mean(stats4.cprom)];
meanMP=[meanMP mean(stats4.maxpr)];
meanRadon=[meanRadon mean(R(:))];
stdRadon=[stdRadon std(R(:))];
% hIntProf=sum(R,1);
% vIntProf=sum(R,2);
% subplot(2,wsplot,3)
% plot(hIntProf)
% subplot(2,wsplot,6);
% plot(vIntProf);
% std(vIntProf);

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
%I = (I - minmaxV(1))*(newMax - newMin)/(minmaxV(2) - minmaxV(1)) + newMin;

%I=V(round(size(V,1)/4):round(3*size(V,1)/4),...
%        round(size(V,1)/4):round(3*size(V,1)/4));
%I = zeros(100,100);
%I(25:75, 25:75) = 1;
figure;
subplot(2,wsplot,1)
imagesc(I),axis image;title('dt08 32 C.mhd');

theta = 0:180;
[R,xp] = radon(I,theta);
[mean(R(:)) std(R(:))];
subplot(2,wsplot,2)
imagesc(theta,xp,R);axis square;
title(strcat(num2str(mean(R(:))),'|',num2str(std(R(:)))));
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
%colormap(hot);

[minmaxV]=[min(min(min(V))) max(max(max(V)))];
[minmaxR]=[min(R(:)) max(R(:))];
%% co-occurrence
[GLCM1 SI] = graycomatrix(I,'GrayLimits',minmaxV,'Offset',occuDir,'NumLevels',nLevels,'Symmetric',true);
stats4=GLCM_Features4(GLCM1,1);
meanAutoc=[meanAutoc mean(stats4.autoc)];
meanCorrp=[meanCorrp mean(stats4.corrp)];
meanContr=[meanContr mean(stats4.contr)];
meanHomog=[meanHomog mean(stats4.homop)];
meanEnerg=[meanEnerg mean(stats4.energ)];
meanEntrop=[meanEntrop mean(stats4.entro)];
meanVar=[meanVar mean(stats4.sosvh)];
meanCS=[meanCS mean(stats4.cshad)];
meanCP=[meanCP mean(stats4.cprom)];
meanMP=[meanMP mean(stats4.maxpr)];
meanRadon=[meanRadon mean(R(:))];
stdRadon=[stdRadon std(R(:))];
% hIntProf=sum(R,1);
% vIntProf=sum(R,2);
% subplot(2,wsplot,3)
% plot(hIntProf)
% subplot(2,wsplot,6);
% plot(vIntProf);
% std(vIntProf);

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
%I = (I - minmaxV(1))*(newMax - newMin)/(minmaxV(2) - minmaxV(1)) + newMin;

%I=V(round(size(V,1)/4):round(3*size(V,1)/4),...
%        round(size(V,1)/4):round(3*size(V,1)/4));
%I = zeros(100,100);
%I(25:75, 25:75) = 1;
figure;
subplot(2,wsplot,1)
imagesc(I),axis image;title('dt08 51 C.mhd');

theta = 0:180;
[R,xp] = radon(I,theta);
[mean(R(:)) std(R(:))];
subplot(2,wsplot,2)
imagesc(theta,xp,R);axis square;
title(strcat(num2str(mean(R(:))),'|',num2str(std(R(:)))));
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
%colormap(hot);

[minmaxV]=[min(min(min(V))) max(max(max(V)))];
[minmaxR]=[min(R(:)) max(R(:))];
%% co-occurrence
[GLCM1 SI] = graycomatrix(I,'GrayLimits',minmaxV,'Offset',occuDir,'NumLevels',nLevels,'Symmetric',true);
stats4=GLCM_Features4(GLCM1,1);
meanAutoc=[meanAutoc mean(stats4.autoc)];
meanCorrp=[meanCorrp mean(stats4.corrp)];
meanContr=[meanContr mean(stats4.contr)];
meanHomog=[meanHomog mean(stats4.homop)];
meanEnerg=[meanEnerg mean(stats4.energ)];
meanEntrop=[meanEntrop mean(stats4.entro)];
meanVar=[meanVar mean(stats4.sosvh)];
meanCS=[meanCS mean(stats4.cshad)];
meanCP=[meanCP mean(stats4.cprom)];
meanMP=[meanMP mean(stats4.maxpr)];
meanRadon=[meanRadon mean(R(:))];
stdRadon=[stdRadon std(R(:))];

% hIntProf=sum(R,1);
% vIntProf=sum(R,2);
% subplot(2,wsplot,3)
% plot(hIntProf)
% subplot(2,wsplot,6);
% plot(vIntProf);
% std(vIntProf);

t = 0:pi/20:2*pi;
R0 = 20; 
x0 = round(100/2); y0 = round(100/2);
xi = R0*cos(t)+x0;
yi = R0*sin(t)+y0;
roimask = poly2mask(xi,yi, 100,100);
pr_r = find(roimask);

subplot(2,wsplot,4)
imshow(roimask)
% [meanVHist,xout3]=hist(R,16);
% bar(xout3, meanVHist, 'b'); 
[R,xp] = radon(roimask,theta);
subplot(2,wsplot,5)
imagesc(theta,xp,R);axis square;
title('R_{\theta} (X\prime)');
xlabel('\theta (degrees)');
ylabel('X\prime');
set(gca,'XTick',0:20:180);
%colormap(hot);

[minmaxR]=[min(R(:)) max(R(:))];
%% co-occurrence
[GLCM1 SI] = graycomatrix(roimask,'GrayLimits',[0 1],'Offset',occuDir,'NumLevels',2,'Symmetric',true);
stats4=GLCM_Features4(GLCM1,1);
meanAutoc=[meanAutoc mean(stats4.autoc)];
meanCorrp=[meanCorrp mean(stats4.corrp)];
meanContr=[meanContr mean(stats4.contr)];
meanHomog=[meanHomog mean(stats4.homop)];
meanEnerg=[meanEnerg mean(stats4.energ)];
meanEntrop=[meanEntrop mean(stats4.entro)];
meanVar=[meanVar mean(stats4.sosvh)];
meanCS=[meanCS mean(stats4.cshad)];
meanCP=[meanCP mean(stats4.cprom)];
meanMP=[meanMP mean(stats4.maxpr)];
meanRadon=[meanRadon mean(R(:))];
stdRadon=[stdRadon std(R(:))];

% hIntProf=sum(R,1);
% vIntProf=sum(R,2);
% subplot(2,wsplot,3)
% plot(hIntProf)
% subplot(2,wsplot,6);
% plot(vIntProf);
% std(vIntProf);


%% Feature Curves

figure
curvesWS=4;
indC=5:7;
indM=4;
indS=1:3;
%% mean Autocorrelation
subplot(curvesWS,curvesWS-1,1)
plot(0:7,meanAutoc),title('mean Autocorrelation');
hold all;
sP=plot(indS-1,meanAutoc(indS),'*c');
set(sP, 'Color',[0 1 0]);
cP=plot(indC-1,meanAutoc(indC),'*r');
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indM-1,meanAutoc(indM),'*m');
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Mean');
hleg1 = legend('Mean','Healthy','Calcified','Mixed');
%hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;

%% mean correlation
subplot(curvesWS,curvesWS-1,2)
plot(0:7,meanCorrp),title('mean correlation');
hold all;
sP=plot(indS-1,meanCorrp(indS),'*c');
set(sP, 'Color',[0 1 0]);
cP=plot(indC-1,meanCorrp(indC),'*r');
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indM-1,meanCorrp(indM),'*m');
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Mean');
hleg1 = legend('Mean','Healthy','Calcified','Mixed');
%hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;

%% mean contrast
subplot(curvesWS,curvesWS-1,3)
plot(0:7,meanContr),title('mean contrast');
hold all;
sP=plot(indS-1,meanContr(indS),'*c');
set(sP, 'Color',[0 1 0]);
cP=plot(indC-1,meanContr(indC),'*r');
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indM-1,meanContr(indM),'*m');
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Mean');
hleg1 = legend('Mean','Healthy','Calcified','Mixed');
%hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;

%% mean homogeneity
subplot(curvesWS,curvesWS-1,4)
plot(0:7,meanHomog),title('mean homogeneity');
hold all;
sP=plot(indS-1,meanHomog(indS),'*c');
set(sP, 'Color',[0 1 0]);
cP=plot(indC-1,meanHomog(indC),'*r');
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indM-1,meanHomog(indM),'*m');
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Mean');
hleg1 = legend('Mean','Healthy','Calcified','Mixed');
%hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;

%% mean energy
subplot(curvesWS,curvesWS-1,5)
plot(0:7,meanEnerg),title('mean energy');
hold all;
sP=plot(indS-1,meanEnerg(indS),'*c');
set(sP, 'Color',[0 1 0]);
cP=plot(indC-1,meanEnerg(indC),'*r');
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indM-1,meanEnerg(indM),'*m');
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Mean');
hleg1 = legend('Mean','Healthy','Calcified','Mixed');
%hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;

%% mean entropy
subplot(curvesWS,curvesWS-1,6)
plot(0:7,meanEntrop),title('mean entropy');
hold all;
sP=plot(indS-1,meanEntrop(indS),'*c');
set(sP, 'Color',[0 1 0]);
cP=plot(indC-1,meanEntrop(indC),'*r');
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indM-1,meanEntrop(indM),'*m');
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Mean');
hleg1 = legend('Mean','Healthy','Calcified','Mixed');
%hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;

%% mean var
subplot(curvesWS,curvesWS-1,7)
plot(0:7,meanVar),title('mean var');
hold all;
sP=plot(indS-1,meanVar(indS),'*c');
set(sP, 'Color',[0 1 0]);
cP=plot(indC-1,meanVar(indC),'*r');
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indM-1,meanVar(indM),'*m');
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Mean');
hleg1 = legend('Mean','Healthy','Calcified','Mixed');
%hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;

%% mean cluster sh
subplot(curvesWS,curvesWS-1,8)
plot(0:7,meanCS),title('mean cluster sh');
hold all;
sP=plot(indS-1,meanCS(indS),'*c');
set(sP, 'Color',[0 1 0]);
cP=plot(indC-1,meanCS(indC),'*r');
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indM-1,meanCS(indM),'*m');
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Mean');
hleg1 = legend('Mean','Healthy','Calcified','Mixed');
%hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;

%% mean cluster prominency
subplot(curvesWS,curvesWS-1,9)
plot(0:7,meanCP),title('mean cluster prominency');
hold all;
sP=plot(indS-1,meanCP(indS),'*c');
set(sP, 'Color',[0 1 0]);
cP=plot(indC-1,meanCP(indC),'*r');
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indM-1,meanCP(indM),'*m');
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Mean');
hleg1 = legend('Mean','Healthy','Calcified','Mixed');
%hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;

%% mean max prob
subplot(curvesWS,curvesWS-1,10)
plot(0:7,meanMP),title('mean max prob');
hold all;
sP=plot(indS-1,meanMP(indS),'*c');
set(sP, 'Color',[0 1 0]);
cP=plot(indC-1,meanMP(indC),'*r');
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indM-1,meanMP(indM),'*m');
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Mean');
hleg1 = legend('Mean','Healthy','Calcified','Mixed');
%hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;

%% mean radon transform
subplot(curvesWS,curvesWS-1,11)
plot(0:7,meanRadon),title('mean RT');
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
subplot(curvesWS,curvesWS-1,12)
plot(0:7,stdRadon),title('std RT');
hold all;
sP=plot(indS-1,stdRadon(indS),'*c');
set(sP, 'Color',[0 1 0]);
cP=plot(indC-1,stdRadon(indC),'*r');
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indM-1,stdRadon(indM),'*m');
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Mean');
hleg1 = legend('Mean','Healthy','Calcified','Mixed');
%hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;
