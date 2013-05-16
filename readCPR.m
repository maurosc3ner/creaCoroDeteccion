%***************************
%*  2D Feature Extractor:
%*  MEAN, STD, Flux
%
%*  Esteban Correa
%* 
%*  2013,3,12
%*  References :
%*  
%*  
%***************************
%
% Features (mean, std, gradMean, gradStd) are extracted with a 2D region of
% interest.
% ROI is defined by the circle equation:
%
% (x-XC)^2  (y-YC)^2 = 1
%******************************


clear all; clc; close all;
datasetDirectory='training/dt08'
addpath src



%% file lecture
refFilename=strcat(datasetDirectory,'/','_creandpuj_reference_rca_path.txt');

cprFilename=strcat(datasetDirectory,'/','_creandpuj_3Dimage_rca.mhd');

reference=load(refFilename);

info = mha_read_header(cprFilename)
V = mha_read_volume(info);

%% Flags
visualDebug=true;

%% roisetup
t = 0:pi/20:2*pi;
radio=4;
R0 = radio/info.PixelDimensions(1); 
x0 = round(size(V,1)/4); y0 = round(size(V,2)/4);
xi = R0*cos(t)+x0;
yi = R0*sin(t)+y0;
roimask = poly2mask(xi,yi, round(size(V,1)/4),round(size(V,2)/4));
pr_r = find(roimask);
meanV=[];
stdV=[];

softM=[];
indSoft=[];
softStd=[];
softVar=[];

calcM=[];
indCalc=[];
calcStd=[];
calcVar=[];

mixM=[];
indMix=[];
mixStd=[];

meanGrad=[];
stdGrad=[];

meaniCf=[];

hx = [-1,0,1];
hy = -hx';


wGrid=2;

if visualDebug
    figure;
end



for slicePos=1:size(V,3)-1
    
    Slice=V(round(size(V,1)/4):round(3*size(V,1)/4),...
        round(size(V,1)/4):round(3*size(V,1)/4),slicePos);
%     sl2Diff=Slice;
%     diff=anisodiff(Slice,20,20,0.25,2);
%     se = strel('disk',round(R0));
%     J1 = imtophat(diff,se);
%     J2 = imtophat(Slice,se);
% 
% %     subplot(2,2,1)
% %     imagesc(Slice), colormap(gray),axis image,title(num2str(reference(slicePos,6)));
% %     subplot(2,2,2)
% %     imagesc(diff), colormap(gray),axis image
% %     subplot(2,2,3)
% %     imagesc(J1), colormap(gray),axis image
% %     subplot(2,2,4)
% %     imagesc(J2), colormap(gray),axis image
% %     
%     Slice=J1;
    
    roimean = mean(Slice(pr_r));
    roistd = std(Slice(pr_r));
    headerImage=strcat('Slice:',num2str(slicePos),...
        ' Mean:',num2str(roimean),...
        ' Std:',num2str(roistd));

    meanV=[meanV, roimean];
    stdV=[stdV, roistd];
    corte='Sano';
    if reference(slicePos,6)==1.0
        %reference(i,6)
        corte='Suave';
        softM=[softM, roimean];
        softStd=[softStd, roistd];
        indSoft=[indSoft, slicePos];
    end
    if reference(slicePos,6)==2.0
        %reference(i,6)
        corte='Calcificada';
        calcM=[calcM, roimean];
        calcStd=[calcStd, roistd];
        indCalc=[indCalc, slicePos];
    end
    if reference(slicePos,6)==3.0
        corte='Mixta';
        mixM=[mixM, roimean];
        mixStd=[mixStd, roistd];
        indMix=[indMix, slicePos];   
    end
    %% Gradient
    %[fx fy]=gradient(Slice);
    fx = imfilter(Slice,hx);
    fy = imfilter(Slice,hy);
    [X Y]=meshgrid(1:size(Slice,1),1:size(Slice,2));
    gradMag=(abs(fx)+abs(fy));
    
    gradroimean = mean(gradMag(pr_r));
    gradroistd = std(gradMag(pr_r));
    
    meanGrad=[meanGrad, gradroimean];
    stdGrad=[stdGrad, gradroistd];
    

    %%Flux
    for i=1:size(Slice,1)
       for j=1:size(Slice,2)
            u=[x0-X(i,j) y0-Y(i,j)];
            u=u./norm(u);
            ux(i,j)=u(1); uy(i,j)=u(2);
            FLUX(i,j)=([fx(i, j) fy(i , j)]*u');
       end
    end
    ind=1;
    for r=1:R0
        k=1;
        for x=x0-r:x0+r
            y= y0+sqrt(   r^2  -   (x-x0)^2    );
            p(k,:)=[x y]; p(k+1,:)=[x y0-sqrt(   r^2  -   (x-x0)^2    )];
            k=k+2;
        end
        Circleflux(ind)=sum(FLUX( sub2ind(size(FLUX), round(p(:,1)),round(p(:,2))) ));
        Circleflux(ind)=Circleflux(ind)/size(p,1);
        %FluxMagnitude(xM,yM)=Circleflux(ind);
        ind=ind+1;
    end
    [iCf idxCf]=max(Circleflux);
    meaniCf=[meaniCf, iCf];
    
    %% Visual Debgging
    if visualDebug
        subplot(wGrid,wGrid,1)
        imagesc(Slice),axis square; colormap(gray), 
        title(strcat('Slice:',corte))
        LineHandler = line(xi,yi,'LineWidth',2,'Color',[.8 0 0]);
        %histogram of original image
        subplot(wGrid,wGrid,2)
        n1=[];
        xout1=[];
        n2=[];
        xout2=[];
        [n1, xout1] =hist(Slice(pr_r),16);
        bar(xout1, n1,'r');grid on;  % Plot in red
        title('Histogram within roi')
        %title(corte);
        %hold on;
        
        subplot(wGrid,wGrid,4)
        imagesc(Slice),axis square; colormap(gray), title(headerImage)
        %subplot(wGrid,wGrid,5)
        [n2, xout2] =hist(gradMag(pr_r),16);
        bar(xout2, n2, 'b');grid on;  % Plot in green
        title('Histogram within gradient magnitude roi')
        %hist(sl2Diff(pr_r))
        %colormap(gray)

        subplot(wGrid,wGrid,3)
        imagesc(Slice), 
        hold on, 
        quiver(X(1:6:end,1:6:end),Y(1:6:end,1:6:end),fx(1:6:end,1:6:end),fy(1:6:end,1:6:end),'Color','b')
        %quiver(X(1:6:end,1:6:end),Y(1:6:end,1:6:end),ux(1:6:end,1:6:end),uy(1:6:end,1:6:end),'Color','r','AutoScale','off')
        hold off,axis square;
        %subplot(wGrid,wGrid,6)
        %[n3, xout3]=hist(gradMag(pr_r))
        %hist(xout3, n3)
        %colormap(gray)
        
%         subplot(wGrid,wGrid,6), imagesc(gradMag), colormap(gray),axis square;
%         subplot(wGrid,wGrid,7),imagesc(FLUX),axis square, title(num2str(iCf));
%         hold on
% 
%         plot([x0 y0], [x0 y0+idxCf],'LineWidth', 2, 'Color','b');
%         hold off
%         % Flux maximun curve
%         subplot(wGrid,wGrid,8)
%          plot(Circleflux), grid on; 
        pause(0.004)
        w = waitforbuttonpress;
    end
    
end


%% Plotting curves
curvesWS=3;
figure, title(datasetDirectory);
subplot(curvesWS,curvesWS-1,1)
plot(meanV),title('Mean Value roi Image');
hold all;
sP=plot(indSoft,softM,'*c')
set(sP, 'Color',[1 0.6471 0]);
cP=plot(indCalc,calcM,'*r')
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indMix,mixM,'*m')
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Mean')
%hleg1 = legend('Mean','Soft','Calcified','Mixed');
hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;

subplot(curvesWS,curvesWS-1,3)
plot(stdV),title('std Value roi Image');
hold all;
sP=plot(indSoft,softStd,'*c')
set(sP, 'Color',[1 0.6471 0]);
cP=plot(indCalc,calcStd,'*r')
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indMix,mixStd,'*m')
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Std')
%hleg1 = legend('Std','Soft','Calcified','Mixed');
hleg1 = legend('Std','Calcified','Mixed');
grid
hold off;

subplot(curvesWS,curvesWS-1,2)
plot(meanGrad),title('Mean Value roi gradient Image');grid;
% hold all;
% sP=plot(indSoft,softM,'*c')
% set(sP, 'Color',[1 0.6471 0]);
% cP=plot(indCalc,calcM,'-*r')
% set(cP, 'Color',[255/255 0/255 0/255]);
% mP=plot(indMix,mixM,'*m')
% set(mP, 'Color',[240/255 128/255 128/255]);
% set(get(gca,'YLabel'),'String','Mean')
% hleg1 = legend('Mean','Soft','Calcified','Mixed');
% grid
% hold off;

subplot(curvesWS,curvesWS-1,4)
plot(stdGrad),title('std Value roi gradient Image');grid;

subplot(curvesWS,curvesWS-1,5)
plot(meaniCf),title('2D Flux');grid
subplot(curvesWS,curvesWS-1,6)
[meanVHist,xout3]=hist(meanV,16);
bar(xout3, meanVHist, 'b'); % Histogram Normalization
title('Histogram along centerline');
%bar(xout3, meanVHist/trapz(xout3,meanVHist), 'b'); % Histogram Normalization
% pdf_g=1/sqrt(2*pi)*exp(-0.5*xout3.^2);%# pdf of the normal distribution
% hold on
% plot(xout3,pdf_g,'r');
% hold off
grid on; 


% hold all;
% sP=plot(indSoft,softStd,'*c')
% set(sP, 'Color',[1 0.6471 0]);
% cP=plot(indCalc,calcStd,'-*r')
% set(cP, 'Color',[255/255 0/255 0/255]);
% mP=plot(indMix,mixStd,'*m')
% set(mP, 'Color',[240/255 128/255 128/255]);
% set(get(gca,'YLabel'),'String','Std')
% hleg1 = legend('Std','Soft','Calcified','Mixed');
% grid
% hold off;


% %%
% windowSize=8;
% lpf = ones(1,windowSize)/windowSize;% variar 40
% %result = conv(meanV, lpf, 'same');
% figure,title('Low pass filter')
% subplot(2,2,1)
% plot(conv(meanV, lpf, 'same')),title('Mean Value roi Image');grid
% subplot(2,2,2)
% plot(conv(stdV, lpf, 'same')),title('std Value roi Image');grid
% subplot(2,2,3)
% plot(conv(meanGrad, lpf, 'same')),title('Mean Value roi gradient Image');grid
% subplot(2,2,4)
% plot(conv(stdGrad, lpf, 'same')),title('std Value roi gradient Image');grid

