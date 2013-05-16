%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Volume interpolatio in CPR Vessels (Circle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: E. Correa, March 18, 2013
% V: CPR.mhd
% Based on:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;
datasetDirectory='training/dt08'
addpath src

%% file lecture
refFilename=strcat(datasetDirectory,'/','_creandpuj_reference_rca_path.txt');

cprFilename=strcat(datasetDirectory,'/','_creandpuj_3Dimage_rca.mhd');

reference=load(refFilename);

info = mha_read_header(cprFilename)
V = mha_read_volume(info);
dims=size(V);

%% finding illness points to interpolate
indInt=find(reference(:,6)~=0);
indGr=[];
for i=1:size(indInt,1)-1
    if(indInt(i)+1~=indInt(i+1))
        disp('diferente grupo')
        indGr=[indGr i]
        indInt(i)
        indInt(i+1)
    end
    
end
indGr=[indGr i+1]

%% cell array to store Vrois
ini=1;
for i=1:size(indGr,2)
    Vroi=V(:,:,ini:indGr(i));
    vector=reshape(Vroi,1,size(Vroi,1)*size(Vroi,2)*size(Vroi,3));
    imgCell{i,1}=vector;
    ini=indGr(i)+1;
end
vectorT=imgCell{1};
vReshaped=reshape(vectorT,dims(1),dims(2),indGr(1));
% checking if they are equal
%sum(vector-vectorT)

%% TOMAR FOTOS!!!!!!!
[x,y,z] = meshgrid(1:1:dims(1), 1:1:dims(2), 1:1:indGr(1));
[xint,yint,zint] = meshgrid(1:dims(1),1:dims(2),linspace(1,indGr(1),indGr(1)+(indGr(1)-1)*3));
vInterpolated =interp3(x,y,z,vReshaped,xint,yint,zint); % vi is x-by-y-by-z+(z-1)*3
colormap(gray)
for slicePos=linspace(1,10,37)
    slice(xint,yint,zint,vInterpolated,[108],[54 108],slicePos), shading flat
    pause(0.5)
end
figure;

circleMask=[1 1 1; 1 1 1; 1 1 1];


%% circle shifting around vessel with square pattern
%coordinates need to be circular!
%xy_1,xy_2...,xy_1
center=[dims(1)/2 dims(2)/2];
step=1.5;
x = [center(1)-step center(1)+step center(1)+step center(1)-step center(1)-step]
y = [center(2)-step center(2)-step center(2)+step center(2)+step center(2)-step]
squareMask = poly2mask(x,y,dims(1),dims(1));

[I,J] = ind2sub(size(squareMask),find(squareMask>0)) 
imagesc(squareMask);
t = 0:pi/20:2*pi;
    radio=4;
    R0 = radio/info.PixelDimensions(1); 
for(i=1:9)
    xi = R0*cos(t)+I(i);
    yi = R0*sin(t)+J(i);
    pause(5);
    LineHandler = line(xi,yi,'LineWidth',1,'Color',[0 1 0]);
end


for slicePos=1:indGr(1)+(indGr(1)-1)*3
   imagesc(vInterpolated(:,:,slicePos)),colormap(jet);

end
   
%% 5 lines circular shifting
A = [ 1:5;6:10; 11:15;16:20;21:25]
%shSeq[ y-1,x-1 y-1,x y-1,x+1;
%        y,x-1   y,x   y,x+1;
%       y+1,x-1 y+1,x y+1,x+1]
shSeq={[-1 -1];[-1 0];[-1 1];...
        [0 -1];[0 0];[0 1];...
        [1 -1];[1 0];[1 1]}; 
for cellPos=1:size(shSeq,1)
    B = circshift(A,shSeq{cellPos})
    cellPos
end
