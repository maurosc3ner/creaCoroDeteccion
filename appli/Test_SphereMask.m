clear all; clc; close all;
datasetDirectory='training/dt07'
addpath src


%% file lecture
refFilename=strcat(datasetDirectory,'/','_creandpuj_reference_rca_path.txt');

cprFilename=strcat(datasetDirectory,'/','_creandpuj_3Dimage_rca.mhd');

reference=load(refFilename);

info = mha_read_header(cprFilename)
V = mha_read_volume(info);

%Sphere mask creation
[dims]=size(V)

XC=round(dims(1)/2);
YC=round(dims(1)/2);
ZC=round(dims(1)/2);
radio=4;     %radius in mm
R0 = radio/info.PixelDimensions(1); 
RX=R0;
RY=R0;
RZ=R0;
[X,Y,Z]=ELLIPSOID(XC,YC,ZC,RX,RY,RZ);

mesh(X,Y,Z)


VM=zeros(dims);
for i=1:dims(1)
    for j=1:dims(2)
        for k=1:dims(3)
            ellipsoidequation=(i-XC)^2/RX^2+...
                (j-YC)^2/RY^2+(k-ZC)^2/RZ^2;
            if ellipsoidequation<1
                %fprintf('Punto por dentro');
                pInside = [pInside sub2ind(dims, i, j, k)];
                VM(i,j,k)=1;
            end
                
           
        end
    end
end

VM=sphereMask(dims,XC,YC,ZC,RX,RY,RZ);
%VM(pInside)=1;

pr_r=find(VM>0);

figure;
for slicePos=1:dims(3)-1
    
    imagesc(VM(:,:,slicePos)),colormap(gray),title(num2str(slicePos));
    pause(0.08);
end
% figure;
% xslice = [dims(1)];
% yslice = dims(2); 
% zslice = [dims(3)];
% slice([1:dims(1)],[1:dims(2)],[1:dims(3)],V,xslice,yslice,zslice)
% colormap jet
%  drawnow
%  
figure;
for slicePos=1:dims(3)-1
    
    Slice=V(:,:,slicePos).*double(VM(:,:,slicePos));
    subplot(2,2,1), imagesc(V(:,:,slicePos)), colormap(gray),axis square;
    subplot(2,2,2),imagesc(Slice),axis square;
    subplot(2,2,3),imagesc(VM(:,:,slicePos)),axis square;
    pause(0.05);
end



pr_r = find(roimask);


[x3,y3,z3]=ELLIPSOID(42,42,42,5,5,5);
mesh(x3,y3,z3)

pInside=[42 42 42];
pOutside=[42 42 42];

if ellipsoidequation>1
    fprintf('Punto por fuera');
elseif ellipsoidequation<1
    fprintf('Punto por dentro');
end



[x,y,z] = meshgrid(-2:.2:2,-2:.25:2,-2:.16:2);
v = x.*exp(-x.^2-y.^2-z.^2);
xslice = [-1.2,.8,2]; yslice = 2; zslice = [-2,0];
slice(x,y,z,v,xslice,yslice,zslice)
colormap hsv


%% Sphere concept Visualization
[xsp,ysp,zsp] = sphere;
slice(x,y,z,v,[-2,0,2],2,-2),title('CPR Vessel')  % Draw some volume boundaries

for i = -3:.2:3
     hsp = surface(xsp+i,ysp,zsp);
     rotate(hsp,[1 0 0],90)
     xd = get(hsp,'XData');
     yd = get(hsp,'YData');
     zd = get(hsp,'ZData');
     delete(hsp)
     hold on
     hslicer = slice(x,y,z,v,xd,yd,zd);
     axis tight ,grid on;
     xlim([-3,3])
     view(-10,35)
     %w = waitforbuttonpress;
     key = input('Press key to continue');
     drawnow
     delete(hslicer)
     hold off
     %pause(0.1);
end



I = checkerboard(8);
PSF = fspecial('gaussian',7,10);
V = .0001;
BlurredNoisy = imnoise(imfilter(I,PSF),'gaussian',0,V);
WT = zeros(size(I));
WT(5:end-4,5:end-4) = 1;
INITPSF = ones(size(PSF));
[J P] = deconvblind(BlurredNoisy,INITPSF,20,10*sqrt(V),WT);
subplot(221);imshow(BlurredNoisy);
title('A = Blurred and Noisy');
subplot(222);imshow(PSF,[]);
title('True PSF');
subplot(223);imshow(J);
title('Deblurred Image');
subplot(224);imshow(P,[]);
title('Recovered PSF');
 