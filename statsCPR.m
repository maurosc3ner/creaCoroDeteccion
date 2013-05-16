clc;clear all; close all;
path='/Users/maurosc3ner/Dropbox/esteban/Creatis/dataset/tr1dt05RCA.mhd' 

img=double(ReadMHD(path));
gx=zeros(size(img));
gy=zeros(size(img));
figure;
for k=1:size(img,3)
    slice=img(:,:,k);
    subplot(2,1,1);
    imagesc(slice);
    [FX,FY]=gradient(slice);
    
    %[FX,FY] = gradient(slice);
    %gx(:,:,k) gy(:,:,k)
    [X,Y]=meshgrid(1:1:60,1:1:60);
    subplot(2,1,2);
    contour(1:1:60,1:1:60,slice),
    hold on,
    quiver(X,Y,FX(1:1:end,1:1:end ),FY(1:1:end,1:1:end ));
    hold off,
    
%     fila=
%     for i=5:34,
%         %change(i)=sum(sum(tore1(:,:,i)-tore1(:,:,i+1))); 
%     end, 
%     %plot(change)
    Fslice=fftshift(fft2(fftshift(slice)));
    fftshow(Fslice);
    pause(0.1);
    %key = input('Press key to continue');
    
end

Rt=27.853703;
Rup_Rsi=1;
for i=1:8
    Rup_Rsi=1+(i/10)
    RSI=((Rup_Rsi+1)*Rt)/Rup_Rsi
    RUP=RSI*Rup_Rsi
end    

%view3dgui(img)
