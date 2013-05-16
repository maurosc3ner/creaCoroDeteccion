%***************************
%*  3D Feature Extractor
%*  Esteban Correa
%* 
%*  2013,3,12
%*  References :
%*  
%*  
%***************************
%
% Features (mean, std, gradMean, gradStd) are extracted with a 3D region of
% interest.
% ROI is defined by the cylinder equation:
%
% (x-XC)^2 + (y-YC)^2 = r^2, z ? [0 L] 
%******************************

clear all; clc; close all;
datasetDirectory='training/dt08'
addpath src
addpath ../ReadData3D_version1k/mha
%% file lecture
refFilename=strcat(datasetDirectory,'/','_creandpuj_reference_rca_path.txt');

cprFilename=strcat(datasetDirectory,'/','_creandpuj_3Dimage_rca.mhd');

reference=load(refFilename);

info = mha_read_header(cprFilename)
V = mha_read_volume(info);

%% Flags
visualDebug=true;

%% Sphere mask creation
[dims]=size(V);



t = pi/4:pi/4:2*pi; 
L=5;
x0 = round(dims(1)/2); y0 = round(dims(2)/2);
%for steps=0:2:20-round(L/2)
for steps=0+20:2:100
    z0=round(L/2)+steps;
    imagesc(V(:,:,z0)),axis square;
    for radio=1:2:5
        %radio=4;
        R0 = radio/info.PixelDimensions(1); 
        yi = R0*cos(t);
        xi = R0*sin(t);
        
        %Mask initialization
        VM=uint8(zeros(dims(1),dims(2),20));
        Slice=uint8(zeros(dims(1),dims(2)));
        
        
        % 2D Sampling pattern
        pr_r=sub2ind(size(Slice),round(x0+xi),round(y0+yi));
        Slice(pr_r)=1;
        % 3D Sampling pattern
        for k=z0-round(L/2)+1:z0+round(L/2)-1;
            VM(:,:,k)=Slice(:,:);
        end
        
        %radial direction points
        ui=[[(x0+xi)-x0]' [(y0+yi)-y0]'];
        %ui=[ col-(col+ux) row-(row+uy)];
        ui=ui./norm(ui);
        ui=ui.*2;
        
        ti=[-ui(:,2) ui(:,1)];
        
        %imagesc(Slice),axis square;
        hold all
        
        %hold on
        plot(x0+xi',y0+yi','d','Color',[1 1 1]),axis square;
        quiver(x0+xi',y0+yi',ui(:,1),ui(:,2))
        quiver(x0+xi',y0+yi',ti(:,1),ti(:,2),'Color','r');
        colormap gray
        hold off
        pause(0.5)
        key=input('key');
    %% Slice visualization
%         [x,y,z] = meshgrid(1:1:dims(1), 1:1:dims(2), 1:1:20);
%         colormap(hot)
%         for slicePos=1:20
%             slice(x,y,z,double(VM),[108],[108],slicePos),shading flat;alpha(0.5);
%             pause(0.2)
%         end
    %% CPR Visualization
%         for slicePos=1:20
%            imagesc(VM(:,:,slicePos)),axis square,colormap gray,title(num2str(slicePos)); 
%            hold on
%            LineHandler = line(xi+x0,yi+y0,'LineWidth',0.5,'Color',[.8 0 0]);
%            hold off;
%            pause(0.5)
%         end


        % for slicePos=1:20
        %    imagesc(VM(:,:,slicePos)),axis square,colormap gray,title(num2str(slicePos)); 
        %    pause(0.5)
        % end
    %     for slicePos=1:20
    %         slice(x,y,z,double(VM),[108],[108],slicePos),shading flat;alpha(0.5);
    %         pause(0.1)
    %     end

    end
    VM=circshift(VM,[0 0 steps]);
end

