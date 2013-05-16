clear all; clc; close all;
datasetDirectory='training/dt01'
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
[dims]=size(V);
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


wGrid=4;

if visualDebug
    figure;
end



for slicePos=1:size(V,3)-1
    
    Slice=V(round(size(V,1)/4):round(3*size(V,1)/4),...
        round(size(V,1)/4):round(3*size(V,1)/4),slicePos);
    


    
    %% Visual Debgging
    if visualDebug
        %Original image
        subplot(2,wGrid,1)
        imagesc(Slice),axis square; colormap(gray), title('original')
        LineHandler = line(xi,yi,'LineWidth',2,'Color',[.8 0 0]);
        %histogram of original image
        subplot(2,wGrid,wGrid+1)
%         n1=[];
%         xout1=[];
%         n2=[];
%         xout2=[];
%         [n1, xout1] =hist(Slice(pr_r),32);
%         bar(xout1, n1,'r');grid on;  % Plot in red
%         title('Histogram within roi')
        diff=anisodiff(Slice,10,20,0.25,2);
        imagesc(diff),axis square;
        %title(corte);
        %hold on;
        
%         subplot(wGrid,wGrid,4)
%         imagesc(Slice),axis square; colormap(gray);
%         %subplot(wGrid,wGrid,5)
%         [n2, xout2] =hist(gradMag(pr_r),16);
%         bar(xout2, n2, 'b');grid on;  % Plot in green
%         title('Histogram within gradient magnitude roi')

        subplot(2,wGrid,2)
        %% Test3 blind deconv testing edge weights
%         WEIGHT = edge(Slice,'sobel');
%         se = strel('disk',2);
%         WEIGHT = 1-double(imdilate(WEIGHT,se));
%         WEIGHT([1:3 end-[0:2]],:) = 0;
%         WEIGHT(:,[1:3 end-[0:2]]) = 0;
%         imshow(WEIGHT);title('Weight array');
%         PSF = fspecial('gaussian');
%         INITPSF = ones(size(PSF));
%         subplot(2,wGrid,3)
%         [J P] = deconvblind(Slice,INITPSF,20);
%         imagesc(J);title('Deblurred Image');
        %% Test2 blind deconv with different PSF
%         PSF = fspecial('gaussian',7,1);
%         UNDERPSF = ones(size(PSF)-4);
%         [J1 P1] = deconvblind(Slice,UNDERPSF);
%         imagesc(J1);title('Deblurring with Undersized PSF'),axis square;
%         subplot(2,wGrid,3)
%         OVERPSF = padarray(UNDERPSF,[4 4],'replicate','both');
%         [J2 P2] = deconvblind(Slice,OVERPSF);
%         imshow(J2);title('Deblurring with Oversized PSF'),axis square, colormap gray;
%         INITPSF = padarray(UNDERPSF,[2 2],'replicate','both');
%         [J3 P3] = deconvblind(Slice,INITPSF);
%         imshow(J3);title('Deblurring with INITPSF'),axis square;
%         
%         subplot(2,wGrid,wGrid+2)
%         imshow(P1,[],'notruesize');,axis square;
%         title('Reconstructed Undersized PSF');
%         subplot(2,wGrid,wGrid+3),axis square;
%         imshow(P2,[],'notruesize');
%         title('Reconstructed Oversized PSF');
%         subplot(2,wGrid,wGrid+3),axis square;
%         imshow(P3,[],'notruesize');
%         title('Reconstructed true PSF');
            
        %% Test1 Wiener
%         interval=[1:size(Slice,1)];
%         interval=interval-round(size(Slice,1)/2);
%         cf = abs(fft2(Slice));
%         surf(interval/size(Slice,1),interval/size(Slice,2),fftshift(cf))
%         shading interp, colormap gray
%         title('abs Slice FFT');
%         
%         %PSNR 40db
%         sigma_u = 10^(-40/20)*abs(1-0);
%         Slice_noise = Slice + sigma_u*randn(size(Slice));
%         h = fspecial('disk',4);
%         hf = fft2(h,size(Slice,1),size(Slice,2));
%         Slice_pinv = real(ifft2((abs(hf) > 0.1).*fft2(Slice_noise)./hf));
%         
%         mse = mean((Slice(:)-Slice_noise(:)).^2)
%         subplot(2,wGrid,4)
        
        
        %imagesc(Slice_pinv),axis square;
%         sigma2_x = var(Slice(:))
%         mean_x = mean(Slice(:))
%         Slice_r = circshift(Slice,[1 0]);
%         Slice_c = circshift(Slice,[0 1]);
%         rho_mat = corrcoef([Slice(:); Slice(:)],[Slice_r(:); Slice_c(:)])
%         rho = rho_mat(1,2);
%         [rr,cc] = ndgrid(interval,interval);
%         r_x = sigma2_x*rho.^sqrt(rr.^2+cc.^2) + mean_x^2;
%         
%         subplot(2,wGrid,wGrid+2)
%         surf(interval,interval,r_x)
%         axis tight
%         shading interp, camlight, colormap jet
%         title('image autocorrelation model approximation')
%         
%         %Wiener deconvolution
%         subplot(2,wGrid,3)
%         h = fspecial('disk',4);
%         %PSNR 40db
%         sigma_u = 10^(-40/10)*abs(1-0);
%         Slice_noise = Slice + sigma_u*randn(size(Slice));
%         Slice_wnr = deconvwnr(Slice_noise,h,numel(Slice)*sigma_u^2./cf);
%         %deconvwnr(Slice,h,sigma_u^2,cf);
%         imagesc(Slice_wnr)
%         colormap(gray)
%         title('restored image using correlation model')
%         mse = mean((Slice(:)-Slice_wnr(:)).^2)
        
        %pause(0.004)
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

