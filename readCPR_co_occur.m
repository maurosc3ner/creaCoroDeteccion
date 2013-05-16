%***************************
%*  2D Feature Extractor:
%*  Co-occurrence features
%  
%*  2013,4,2
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
datasetDirectory='training/dt09'
addpath src
addpath src/GLCM_Features4


%% file lecture
refFilename=strcat(datasetDirectory,'/','_creandpuj_reference_rca_path.txt');

cprFilename=strcat(datasetDirectory,'/','_creandpuj_3Dimage_rca.mhd');

reference=load(refFilename);

info = mha_read_header(cprFilename)
V = mha_read_volume(info);
[dims]=size(V);
[minmaxV]=[min(min(min(V))) max(max(max(V)))]
%% Flags
visualDebug=false;

%% roisetup
t = 0:pi/20:2*pi;
radio=4;
R0 = radio/info.PixelDimensions(1); 
x0 = round(R0)+1; y0 = round(R0)+1;
xi = R0*cos(t)+x0;
yi = R0*sin(t)+y0;
roimask = poly2mask(xi,yi, round(R0)*2+1,round(R0)*2+1);
pr_r = find(roimask);

radio2=3;
R1 = (radio2/info.PixelDimensions(1))/2; 
xi2 = R1*cos(t)+x0;
yi2 = R1*sin(t)+y0;
roimask2 = poly2mask(xi2,yi2, round(dims(1)/4),round(dims(2)/4));
pr_r2 = find(roimask2);

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

meanContr=[];
meanHomog=[];
meanEnerg=[];
meanEntro=[];
wGrid=2;

if visualDebug
    figure;
end
%co-occur 1mm
occOff=round(1/info.PixelDimensions(1));
radio/info.PixelDimensions(1)
occuDir=[0 occOff;   
        -occOff occOff;    
        -occOff 0;    
        -occOff -occOff];    
%         0 -occOff;     
%         occOff -occOff;    
%         occOff 0;     
%         occOff occOff]

for slicePos=1:dims(3)-1
    
    if(mod(dims(1),2)==0)
        Slice=V(round(dims(1)/2-R0):round(dims(1)/2+R0),...
            round(dims(2)/2-R0):round(dims(2)/2+R0),slicePos);
    else
        Slice=V(round(dims(1)/2-R0):round(dims(1)/2+R0)+1,...
            round(dims(2)/2-R0):round(dims(2)/2+R0)+1,slicePos);
    end
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
    [fx fy]=gradient(Slice);
%     fx = imfilter(Slice,hx);
%     fy = imfilter(Slice,hy);
    [X Y]=meshgrid(1:size(Slice,1),1:size(Slice,2));
    %gradMag=(abs(fx)+abs(fy));
    
    [GLCM1 SI] = graycomatrix(Slice,'GrayLimits',minmaxV,'Offset',occuDir,'NumLevels',64,'Symmetric',true);
    %stats = graycoprops(GLCM1,{'contrast','homogeneity','energy',});
%     stats.Contrast
%     stats.Homogeneity
%     stats.Energy
    stats4 = GLCM_Features4(GLCM1,1);
    contr=max(stats4.contr);
    homog=max(stats4.homop);
    entrop=max(stats4.entro);
    %entro=entropy(Slice(pr_r))
    
    meanContr=[meanContr contr];
    meanHomog=[meanHomog homog];
    meanEnerg=[meanEnerg entrop];
    key=input('input');
    %% Visual Debgging
    if visualDebug
        subplot(wGrid,wGrid,1)
        imagesc(Slice),axis square; colormap(gray), title('original')
        LineHandler = line(xi,yi,'LineWidth',2,'Color',[.8 0 0]);
        %histogram of original image
        subplot(wGrid,wGrid,2)

        %title(corte);
        %hold on;
        
        subplot(wGrid,wGrid,4)
        imagesc(Slice),axis square; colormap(gray), title(num2str(contr))
        %subplot(wGrid,wGrid,5)

        %hist(sl2Diff(pr_r))
        %colormap(gray)

        subplot(wGrid,wGrid,3)

        %subplot(wGrid,wGrid,6)
        %[n3, xout3]=hist(gradMag(pr_r))
        %hist(xout3, n3)
        %colormap(gray)
        

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

%% mean Contrast
subplot(curvesWS,curvesWS-1,2)
plot(meanContr),title('Co-Occurrence mean contrast');
hold all;
sP=plot(indSoft,meanContr(indSoft),'*c')
set(sP, 'Color',[1 0.6471 0]);
cP=plot(indCalc,meanContr(indCalc),'*r')
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indMix,meanContr(indMix),'*m')
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Mean')
%hleg1 = legend('Mean','Soft','Calcified','Mixed');
hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;

%% mean Homogeneity
subplot(curvesWS,curvesWS-1,4)
plot(meanHomog),title('Co-Occurrence mean homogeneity');
hold all;
sP=plot(indSoft,meanHomog(indSoft),'*c')
set(sP, 'Color',[1 0.6471 0]);
cP=plot(indCalc,meanHomog(indCalc),'*r')
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indMix,meanHomog(indMix),'*m')
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Mean')
%hleg1 = legend('Mean','Soft','Calcified','Mixed');
hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;

%% mean Homogeneity
subplot(curvesWS,curvesWS-1,6)
plot(meanEnerg),title('Co-Occurrence mean energy');
hold all;
sP=plot(indSoft,meanEnerg(indSoft),'*c')
set(sP, 'Color',[1 0.6471 0]);
cP=plot(indCalc,meanEnerg(indCalc),'*r')
set(cP, 'Color',[255/255 0/255 0/255]);
mP=plot(indMix,meanEnerg(indMix),'*m')
set(mP, 'Color',[240/255 128/255 128/255]);
set(get(gca,'YLabel'),'String','Mean')
%hleg1 = legend('Mean','Soft','Calcified','Mixed');
hleg1 = legend('Mean','Calcified','Mixed');
grid
hold off;
