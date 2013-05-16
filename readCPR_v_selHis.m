%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   ROI Histogram analysis in CPR Vessels (Circle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: E. Correa, March 17, 2013
% V: CPR.mhd
% Based on:
% Teßmann M., Vega-Higuera F., Bischoff B., Hausleiter J., Greiner G.,
% Automatic detection and quantification of coronary calcium on 3D CT 
% angiography data, August 2010.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;
datasetDirectory='training/dt11'
addpath src

%% file lecture
refFilename=strcat(datasetDirectory,'/','_creandpuj_reference_rca_path.txt');

cprFilename=strcat(datasetDirectory,'/','_creandpuj_3Dimage_rca.mhd');

reference=load(refFilename);

info = mha_read_header(cprFilename)
V = mha_read_volume(info);
dims=size(V);
%HU=GV-1024
%V=V-1024;
%% Flags
visualDebug=false;

% %% roisetup
t = 0:pi/20:2*pi;
radio=4;
R0 = radio/info.PixelDimensions(1); 
x0 = round(dims(1)/2); y0 = round(dims(2)/2);
xi = R0*cos(t)+x0;
yi = R0*sin(t)+y0;
roimask = poly2mask(xi,yi, round(dims(1)),round(dims(2)));
pr_r = find(roimask);

Vroi=zeros(size(pr_r),dims(3));

%% Pre-processing TeBmann10
centers=linspace(-1024,1586,128);
[volumeHist xout4]=hist(V(:),128);
sum(volumeHist)
delta_xout4 = xout4(2) - xout4(1);
%histogram normalization by AREA
volumeHist=volumeHist./sum(volumeHist(:)*delta_xout4);
bar(xout4, volumeHist, 'b')
%hist(V(:),xout4);
title('Volume Histogram');
grid on; 

%% smooth by Gaussian 1D
sigma = 0.5;
gsize = 9;
x = linspace(-gsize / 2, gsize / 2, gsize);        %range
gaussFilter = (1/(sqrt(2*pi)*sigma))*exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normal
vHFilt = conv (volumeHist', gaussFilter, 'same');
hold all
plot(xout4,vHFilt, 'LineWidth', 2, 'Color','g')
hold off

%% Threshold search
lastMax=-10000000;
idlm=1;
maxValue=1586;
minDistance=2610;
[pks,locs] = findpeaks(vHFilt);
%re-findpeaks with minpeakheight > mean(heights)
[pks,locs] = findpeaks(vHFilt,'MINPEAKHEIGHT',mean(pks));
hold all
plot(xout4(locs),pks+0.00001,'k^','markerfacecolor',[1 0 0])
hold off

vHFiltGRAD=gradient(vHFilt(locs(end):end));

finalInterval=find(vHFiltGRAD<-0.0001);
finalInterval=finalInterval+locs(end);
hold all
plot(xout4(finalInterval),vHFilt(finalInterval)+0.00001,'kv','markerfacecolor',[1 1 0])
hold off
plaqueTh=[xout4(finalInterval(1)) xout4(finalInterval(end))]

calcM=[];
indCalc=[];
mycalcM=[];
myindCalc=[];

wGrid=2;

if visualDebug
    figure;
end
% 
% roiSlice=round(size(V,1)/4):round(3*size(V,1)/4);
% 

for slicePos=1:dims(3)-1
    
    Slice=V(:,:,slicePos);
    % rule 1
    %idBW=find(Slice(:)>plaqueTh(1) & Slice(:)>plaqueTh(2));
    % rule 2
    idBW=find(Slice(:)>plaqueTh(2));
    BW=zeros(size(Slice));
    BW(idBW)=1;
    BW=BW.*roimask;
    imResult=BW.*Slice;
    area=sum(imResult(:))/(size(idBW,1)+0.0000001);
    mycalcM=[mycalcM, area];
    myindCalc=[myindCalc, slicePos];
    corte='Sano';
    if reference(slicePos,6)==1.0
        corte='Suave';
    elseif reference(slicePos,6)==2.0
        corte='Calcificada';
        indCalc=[indCalc, slicePos];
    elseif reference(slicePos,6)==3.0
        corte='Mixta';
    end

    % Visual Debgging
    if visualDebug
        subplot(wGrid,wGrid,1)
        %figure(1)
        imagesc(Slice),axis square; colormap(gray), title(strcat('original :',num2str(slicePos)));
        LineHandler = line(xi,yi,'LineWidth',2,'Color',[.8 0 0]);
        %histogram of original image
        subplot(wGrid,wGrid,2)
        imagesc(BW),axis square; colormap(gray), 
        title(strcat('threshold Area:',num2str(area)));
        subplot(wGrid,wGrid,3)
        imagesc(imResult),axis square; colormap(gray), title(strcat('threshold ',corte));
        
%         imcontrast()
%         subplot(wGrid,wGrid,2)
%         n1=[];
%         xout1=[];
%         n2=[];
%         xout2=[];
%         [n1, xout1] =hist(Slice(pr_r),16);
%         bar(xout1, n1,'r');grid on;  % Plot in red
%         title('Histogram within roi')
%         %title(corte);
%         %hold on;
%         
%         subplot(wGrid,wGrid,4)
%         imagesc(Slice),axis square; colormap(gray), title(headerImage)
%         %subplot(wGrid,wGrid,5)
%         [n2, xout2] =hist(gradMag(pr_r),16);
%         bar(xout2, n2, 'b');grid on;  % Plot in green
%         title('Histogram within gradient magnitude roi')
%         %hist(sl2Diff(pr_r))
%         %colormap(gray)
% 
%         subplot(wGrid,wGrid,3)
%         imagesc(Slice), 
%         hold on, 
%         quiver(X(1:6:end,1:6:end),Y(1:6:end,1:6:end),fx(1:6:end,1:6:end),fy(1:6:end,1:6:end),'Color','b')
%         %quiver(X(1:6:end,1:6:end),Y(1:6:end,1:6:end),ux(1:6:end,1:6:end),uy(1:6:end,1:6:end),'Color','r','AutoScale','off')
%         hold off,axis square;
%         subplot(wGrid,wGrid,6)
%         [n3, xout3]=hist(gradMag(pr_r))
%         hist(xout3, n3)
%         colormap(gray)
        
        pause(0.004)
        key=input('input');
        %w = waitforbuttonpress;
    end
    
end


%% Plotting curves
curvesWS=2;
figure, title(datasetDirectory);
subplot(curvesWS,curvesWS-1,1)
plot(mycalcM),title('Area Value roi Image');
hold all;
% sP=plot(myindCalc,mycalcM,'*c')
% set(sP, 'Color',[1 0.6471 0]);
cP=plot(indCalc,mycalcM(indCalc),'*r')
set(cP, 'Color',[255/255 0/255 0/255]);
set(get(gca,'YLabel'),'String','area')
%hleg1 = legend('Mean','Soft','Calcified','Mixed');
hleg1 = legend('area','Calcified','Mixed');
grid
hold off;

% curvesWS=3;
% figure, title(datasetDirectory);
% subplot(curvesWS,curvesWS-1,1)
% plot(meanV),title('Mean Value roi Image');
% hold all;
% sP=plot(indSoft,softM,'*c')
% set(sP, 'Color',[1 0.6471 0]);
% cP=plot(indCalc,calcM,'*r')
% set(cP, 'Color',[255/255 0/255 0/255]);
% mP=plot(indMix,mixM,'*m')
% set(mP, 'Color',[240/255 128/255 128/255]);
% set(get(gca,'YLabel'),'String','Mean')
% %hleg1 = legend('Mean','Soft','Calcified','Mixed');
% hleg1 = legend('Mean','Calcified','Mixed');
% grid
% hold off;
% 
% subplot(curvesWS,curvesWS-1,3)
% plot(stdV),title('std Value roi Image');
% hold all;
% sP=plot(indSoft,softStd,'*c')
% set(sP, 'Color',[1 0.6471 0]);
% cP=plot(indCalc,calcStd,'*r')
% set(cP, 'Color',[255/255 0/255 0/255]);
% mP=plot(indMix,mixStd,'*m')
% set(mP, 'Color',[240/255 128/255 128/255]);
% set(get(gca,'YLabel'),'String','Std')
% %hleg1 = legend('Std','Soft','Calcified','Mixed');
% hleg1 = legend('Std','Calcified','Mixed');
% grid
% hold off;
% 
% subplot(curvesWS,curvesWS-1,2)
% plot(meanGrad),title('Mean Value roi gradient Image');grid;
% % hold all;
% % sP=plot(indSoft,softM,'*c')
% % set(sP, 'Color',[1 0.6471 0]);
% % cP=plot(indCalc,calcM,'-*r')
% % set(cP, 'Color',[255/255 0/255 0/255]);
% % mP=plot(indMix,mixM,'*m')
% % set(mP, 'Color',[240/255 128/255 128/255]);
% % set(get(gca,'YLabel'),'String','Mean')
% % hleg1 = legend('Mean','Soft','Calcified','Mixed');
% % grid
% % hold off;
% 
% subplot(curvesWS,curvesWS-1,4)
% plot(stdGrad),title('std Value roi gradient Image');grid;
% 
% 
% subplot(curvesWS,curvesWS-1,5)
% plot(meaniCf),title('2D Flux');grid
% 
% figure;
% %subplot(curvesWS,curvesWS-1,6)
% [meanVHist,xout3]=hist(meanV,32);
% bar(xout3, meanVHist, 'b'); % Histogram Normalization
% title('Histogram along centerline');
% 
% %distance between center bin
% w = xout3(2)-xout3(1);
% %interval creation
% t = linspace(xout3(1)-w/2,xout3(end)+w/2,length(xout3)+1);
% dt = diff(t);
% Fvals = cumsum([0,meanVHist.*dt]);
% %curve
% F = spline(t, [0, Fvals, 0]);
% 
% DF = fnder(F);  % computes its first derivative
% % set(h,'String','h(i)')
% % set(tL,'String','t(i)')
% % set(tR,'String','t(i+1)')
% hold on
% fnplt(DF, 'r', 2)
% hold off
% ylims = ylim; ylim([0,ylims(2)]);
% 
% %bar(xout3, meanVHist/trapz(xout3,meanVHist), 'b'); % Histogram Normalization
% % pdf_g=1/sqrt(2*pi)*exp(-0.5*xout3.^2);%# pdf of the normal distribution
% % hold on
% % plot(xout3,pdf_g,'r');
% % hold off
% grid on; 


% figure;
% %subplot(curvesWS,curvesWS-1,6)


